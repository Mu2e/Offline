// Andrei Gaponenko, 2015

#include <string>
#include <memory>
#include <array>
#include <vector>
#include <algorithm>

#include "cetlib/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Utilities/InputTag.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "G4Track.hh"
#include "G4Step.hh"

#include "Mu2eG4/inc/IMu2eG4Cut.hh"
#include "Mu2eG4/inc/Mu2eG4ResourceLimits.hh"
#include "MCDataProducts/inc/ProcessCode.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
//#include "Mu2eG4/inc/getPhysicalVolumeOrThrow.hh"


namespace mu2e {
  namespace Mu2eG4Cuts {
    using namespace std;
    typedef std::vector<fhicl::ParameterSet> PSVector;

    //================================================================
    // A common implementation for some of the required IMu2eG4Cut methods
    class IOHelper: virtual public IMu2eG4Cut {
    public:
      virtual void declareProducts(art::EDProducer *parent) override;
      virtual void finishConstruction(const CLHEP::Hep3Vector& mu2eOriginInWorld) override;
      virtual void beginEvent(const art::Event& evt, const SimParticleHelper& spHelper) override;
      virtual void put(art::Event& event) override;

    protected:
      explicit IOHelper(const fhicl::ParameterSet& pset, const Mu2eG4ResourceLimits& mu2elimits)
        : outputName_(pset.get<string>("write", ""))
        , preSimulatedHitTag_(pset.get<art::InputTag>("preSimulatedHits", art::InputTag()))
        , spHelper_()
        , mu2elimits_(&mu2elimits)
        , overflowWarningPrinted_(false)
      {
        if((preSimulatedHitTag_ != art::InputTag()) && outputName_.empty()) {
          throw cet::exception("CONFIG")<<"Mu2eG4 Cuts: preSimulatedHits = "<<preSimulatedHitTag_
                                        <<" has no effect as no output is specified via 'write'\n";
        }
      }

      std::string outputName_;
      std::unique_ptr<StepPointMCCollection> output_;
      art::InputTag preSimulatedHitTag_;
      CLHEP::Hep3Vector mu2eOrigin_;
      const SimParticleHelper *spHelper_;
      const Mu2eG4ResourceLimits *mu2elimits_;
      bool overflowWarningPrinted_; // in the current event

      void addHit(const G4Step *aStep, ProcessCode endCode);
    };

    void IOHelper::declareProducts(art::EDProducer *parent) {
      if(!outputName_.empty()) {
        parent->produces<StepPointMCCollection>(outputName_);
      }
    }

    void IOHelper::beginEvent(const art::Event& evt, const SimParticleHelper& spHelper) {
      spHelper_ = &spHelper;
      overflowWarningPrinted_ = false;
      if(!outputName_.empty()) {
        output_ = make_unique<StepPointMCCollection>();
        std::cout<<"IOHelper: Creating output collection "<<outputName_
                 <<" ptr = "<<output_.get()
                 <<std::endl;

        if(preSimulatedHitTag_ != art::InputTag()) {
          const auto& inhits = evt.getValidHandle<StepPointMCCollection>(preSimulatedHitTag_);
          output_->reserve(inhits->size());
          for(const auto& hit: *inhits) {
            output_->emplace_back(hit);
          }
        }

      }
    }

    void IOHelper::put(art::Event& evt) {
      if(output_) {
        std::cout<<"IOHelper: recording output collection "<<outputName_<<std::endl;
        evt.put(std::move(output_), outputName_);
      }
    }

    void IOHelper::finishConstruction(const CLHEP::Hep3Vector& mu2eOriginInWorld) {
      mu2eOrigin_ = mu2eOriginInWorld;
    }

    void IOHelper::addHit(const G4Step *aStep, ProcessCode endCode) {
      if(output_->size() < mu2elimits_->maxStepPointCollectionSize()) {

        // The point's coordinates are saved in the mu2e coordinate system.
        output_->
          push_back(StepPointMC(spHelper_->particlePtr(aStep->GetTrack()),
                                aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(),
                                aStep->GetTotalEnergyDeposit(),
                                aStep->GetNonIonizingEnergyDeposit(),
                                aStep->GetPreStepPoint()->GetGlobalTime(),
                                aStep->GetPreStepPoint()->GetProperTime(),
                                aStep->GetPreStepPoint()->GetPosition() - mu2eOrigin_,
                                aStep->GetPreStepPoint()->GetMomentum(),
                                aStep->GetStepLength(),
                                endCode
                                ));
      }
      else {
        if(!overflowWarningPrinted_) {
          overflowWarningPrinted_ = true;
          mf::LogWarning("G4") << "Maximum number of entries reached in output collection "
                               << outputName_ << ": " << output_->size() << endl;
        }
      }
    }


    //================================================================
    class Union: virtual public IMu2eG4Cut,
                 public IOHelper
    {
    public:
      virtual bool steppingActionCut(const G4Step  *step);
      virtual bool stackingActionCut(const G4Track *trk);

      // Sequences need a different implementation
      virtual void declareProducts(art::EDProducer *parent) override;
      virtual void beginEvent(const art::Event& evt, const SimParticleHelper& spHelper) override;
      virtual void put(art::Event& event) override;
      virtual void finishConstruction(const CLHEP::Hep3Vector& mu2eOriginInWorld) override;

      explicit Union(const fhicl::ParameterSet& pset, const Mu2eG4ResourceLimits& lim);
    private:
      std::vector<std::unique_ptr<IMu2eG4Cut> > cuts_;
    };

    Union::Union(const fhicl::ParameterSet& pset, const Mu2eG4ResourceLimits& lim)
      : IOHelper(pset, lim)
    {
      PSVector pars = pset.get<PSVector>("pars");
      for(const auto& p: pars) {
        cuts_.emplace_back(createMu2eG4Cuts(p, lim));
      }
    }

    bool Union::steppingActionCut(const G4Step *step) {
      bool result = false;
      for(const auto& cut : cuts_) {
        if(cut->steppingActionCut(step)) {
          result = true;
          if(output_) {
            addHit(step, ProcessCode::mu2eKillerVolume);
          }
          break;
        }
      }
      return result;
    }

    bool Union::stackingActionCut(const G4Track *trk) {
      bool result = false;
      for(const auto& cut : cuts_) {
        if(cut->stackingActionCut(trk)) {
          result = true;
          break;
        }
      }
      return result;
    }

    void Union::declareProducts(art::EDProducer *parent) {
      IOHelper::declareProducts(parent);
      for(auto& cut: cuts_) {
        cut->declareProducts(parent);
      }
    }

    void Union::finishConstruction(const CLHEP::Hep3Vector& mu2eOriginInWorld) {
      IOHelper::finishConstruction(mu2eOriginInWorld);
      for(auto& cut: cuts_) {
        cut->finishConstruction(mu2eOriginInWorld);
      }
    }

    void Union::beginEvent(const art::Event& evt, const SimParticleHelper& spHelper) {
      IOHelper::beginEvent(evt, spHelper);
      for(auto& cut: cuts_) {
        cut->beginEvent(evt, spHelper);
      }
    }

    void Union::put(art::Event& evt) {
      IOHelper::put(evt);
      for(auto& cut: cuts_) {
        cut->put(evt);
      }
    }

    //================================================================
    class Intersection: virtual public IMu2eG4Cut,
                        public IOHelper
    {
    public:
      virtual bool steppingActionCut(const G4Step  *step);
      virtual bool stackingActionCut(const G4Track *trk);

      // Sequences need a different implementation
      virtual void declareProducts(art::EDProducer *parent) override;
      virtual void beginEvent(const art::Event& evt, const SimParticleHelper& spHelper) override;
      virtual void put(art::Event& event) override;
      virtual void finishConstruction(const CLHEP::Hep3Vector& mu2eOriginInWorld) override;

      explicit Intersection(const fhicl::ParameterSet& pset, const Mu2eG4ResourceLimits& lim);
    private:
      std::vector<std::unique_ptr<IMu2eG4Cut> > cuts_;
    };

    Intersection::Intersection(const fhicl::ParameterSet& pset, const Mu2eG4ResourceLimits& lim)
      : IOHelper(pset, lim)
    {
      PSVector pars = pset.get<PSVector>("pars");
      for(const auto& p: pars) {
        cuts_.emplace_back(createMu2eG4Cuts(p, lim));
      }
    }

    bool Intersection::steppingActionCut(const G4Step *step) {
      bool result = true;
      for(const auto& cut : cuts_) {
        if(!cut->steppingActionCut(step)) {
          result = false;
          break;
        }
      }
      if(result && output_) {
        addHit(step, ProcessCode::mu2eKillerVolume);
      }

      return result;
    }

    bool Intersection::stackingActionCut(const G4Track *trk) {
      bool result = true;
      for(const auto& cut : cuts_) {
        if(!cut->stackingActionCut(trk)) {
          result = false;
          break;
        }
      }
      return result;
    }

    void Intersection::declareProducts(art::EDProducer *parent) {
      IOHelper::declareProducts(parent);
      for(auto& cut: cuts_) {
        cut->declareProducts(parent);
      }
    }

    void Intersection::finishConstruction(const CLHEP::Hep3Vector& mu2eOriginInWorld) {
      IOHelper::finishConstruction(mu2eOriginInWorld);
      for(auto& cut: cuts_) {
        cut->finishConstruction(mu2eOriginInWorld);
      }
    }

    void Intersection::beginEvent(const art::Event& evt, const SimParticleHelper& spHelper) {
      IOHelper::beginEvent(evt, spHelper);
      for(auto& cut: cuts_) {
        cut->beginEvent(evt, spHelper);
      }
    }

    void Intersection::put(art::Event& evt) {
      IOHelper::put(evt);
      for(auto& cut: cuts_) {
        cut->put(evt);
      }
    }

    //================================================================
    class Plane: virtual public IMu2eG4Cut,
                 public IOHelper
    {
    public:
      virtual bool steppingActionCut(const G4Step  *step);
      virtual bool stackingActionCut(const G4Track *trk);

      explicit Plane(const fhicl::ParameterSet& pset, const Mu2eG4ResourceLimits& lim);
    private:
      std::array<double,3> normal_;
      double offset_;

      bool cut_impl(const CLHEP::Hep3Vector& pos);
    };

    Plane::Plane(const fhicl::ParameterSet& pset, const Mu2eG4ResourceLimits& lim)
      : IOHelper(pset, lim)
      , offset_()
    {
      // FIXME: use pset.get<array>() when it is available
      vector<double> n{pset.get<vector<double> >("normal")};
      if(n.size() != 3) {
        throw std::runtime_error("SteppingCut::Plane(): normal should be a vector of 3 doubles. Error in pset = "+pset.to_string());
      }
      //std::copy(n.begin(), n.end(), &normal_[0]);
      std::copy(n.begin(), n.end(), normal_.begin());

      vector<double> x0{pset.get<vector<double> >("point")};
      if(x0.size() != 3) {
        throw cet::exception("CONFIG")<<"SteppingCut::Plane(): normal should be a vector of 3 doubles. "
                                      <<"Error in pset = "<<pset.to_string()<<"\n";

      }
      // the cut:
      //             (x-x0)*normal >= 0
      // rewrite as
      //
      //              x*normal >= x0*normal =: offset

      offset_ = std::inner_product(normal_.begin(), normal_.end(), x0.begin(), 0.);
    }

    bool Plane::cut_impl(const CLHEP::Hep3Vector& pos) {
      const bool result =
        (pos.x()*normal_[0] + pos.y()*normal_[1] + pos.z()*normal_[2] >= offset_);

      return result;
    }

    bool Plane::steppingActionCut(const G4Step *step) {
      const CLHEP::Hep3Vector& pos = step->GetPostStepPoint()->GetPosition();
      const bool result = cut_impl(pos);
      if(result && output_) {
        addHit(step, ProcessCode::mu2eKillerVolume);
      }
      return result;
    }

    bool Plane::stackingActionCut(const G4Track *trk) {
      const CLHEP::Hep3Vector& pos = trk->GetPosition();
      return cut_impl(pos);
    }

    //================================================================
    class Constant: virtual public IMu2eG4Cut,
                    public IOHelper
    {
    public:
      virtual bool steppingActionCut(const G4Step  *step);
      virtual bool stackingActionCut(const G4Track *trk);

      explicit Constant(bool val, const Mu2eG4ResourceLimits& lim);
      explicit Constant(const fhicl::ParameterSet& pset, const Mu2eG4ResourceLimits& lim);
    private:
      bool value_;
    };

    Constant::Constant(const fhicl::ParameterSet& pset, const Mu2eG4ResourceLimits& lim)
      : IOHelper(pset, lim)
      , value_{pset.get<double>("value")}
    {}

    Constant::Constant(bool val, const Mu2eG4ResourceLimits& lim) : IOHelper(fhicl::ParameterSet(), lim), value_(val) {}

    bool Constant::steppingActionCut(const G4Step *step) {
      if(output_) {
        addHit(step, ProcessCode::mu2eKillerVolume);
      }
      return value_;
    }

    bool Constant::stackingActionCut(const G4Track *) {
      return value_;
    }

    //================================================================
    // FIXME: more cuts:
    //    VolumeCut, PDGIdCut, PrimaryOnly
    //     if ( _primaryOnly ){
    //       if ( trk->GetParentID() != 0 ) {
    //         return fKill;
    //       }
    //     }

    //================================================================
  } // end namespace Mu2eG4Cuts

  //================================================================
  std::unique_ptr<IMu2eG4Cut> createMu2eG4Cuts(const fhicl::ParameterSet& pset, const Mu2eG4ResourceLimits& lim) {
    using namespace Mu2eG4Cuts;

    if(pset.is_empty()) return make_unique<Constant>(false, lim); // no cuts

    const string cuttype =  pset.get<string>("type");

    if(cuttype == "union") return make_unique<Union>(pset, lim);
    if(cuttype == "intersection") return make_unique<Intersection>(pset, lim);
    if(cuttype == "plane") return make_unique<Plane>(pset, lim);
    if(cuttype == "const") return make_unique<Constant>(pset, lim);

    throw cet::exception("CONFIG")<< "mu2e::createMu2eG4Cuts(): can not parse pset = "<<pset.to_string()<<"\n";
  }

} // end namespace mu2e
