// Andrei Gaponenko, 2015

#include <string>
#include <memory>
#include <array>
#include <vector>
#include <algorithm>
#include <numeric>

#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ProducesCollector.h"
#include "art/Framework/Core/ConsumesCollector.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "Geant4/G4Track.hh"
#include "Geant4/G4Step.hh"
#include "Geant4/G4VProcess.hh"

#include "Mu2eG4/inc/IMu2eG4Cut.hh"
#include "Mu2eG4/inc/Mu2eG4ResourceLimits.hh"
#include "MCDataProducts/inc/ProcessCode.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "Mu2eG4/inc/getPhysicalVolumeOrThrow.hh"
#include "DataProducts/inc/PDGCode.hh"

#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"

namespace mu2e {
  namespace Mu2eG4Cuts {
    using namespace std;
    typedef std::vector<fhicl::ParameterSet> PSVector;

    //================================================================
    // A common implementation for some of the required IMu2eG4Cut methods
    class IOHelper: virtual public IMu2eG4Cut {
    public:
      virtual void declareProducts(art::ProducesCollector& pc, art::ConsumesCollector& cc) override;
      virtual void finishConstruction(const CLHEP::Hep3Vector& mu2eOriginInWorld) override;
      virtual void beginEvent(const art::Event& evt, const SimParticleHelper& spHelper) override;
      virtual void put(art::Event& event) override;
      virtual void deleteCutsData() override;

    protected:
      explicit IOHelper(const fhicl::ParameterSet& pset, const Mu2eG4ResourceLimits& mu2elimits)
        : steppingOutputName_(pset.get<string>("write", ""))
        , preSimulatedHitTag_(pset.get<art::InputTag>("preSimulatedHits", art::InputTag()))
        , spHelper_()
        , mu2elimits_(&mu2elimits)
        , overflowWarningPrinted_(false)
      {}

      std::string steppingOutputName_;
      std::unique_ptr<StepPointMCCollection> steppingOutput_;
      art::InputTag preSimulatedHitTag_;
      CLHEP::Hep3Vector mu2eOrigin_;
      const SimParticleHelper *spHelper_;
      const Mu2eG4ResourceLimits *mu2elimits_;
      bool overflowWarningPrinted_; // in the current event

      void addHit(const G4Step *aStep);
    };

    void IOHelper::declareProducts(art::ProducesCollector& pc, art::ConsumesCollector& cc) {

      if(!steppingOutputName_.empty()) {
        if(preSimulatedHitTag_ != art::InputTag()) {
          cc.consumes<StepPointMCCollection>(preSimulatedHitTag_);
        }

       pc.produces<StepPointMCCollection>(steppingOutputName_);
      }
    }

    void IOHelper::beginEvent(const art::Event& evt, const SimParticleHelper& spHelper) {
      spHelper_ = &spHelper;
      overflowWarningPrinted_ = false;
      if(!steppingOutputName_.empty()) {
        steppingOutput_ = make_unique<StepPointMCCollection>();
        if(preSimulatedHitTag_ != art::InputTag()) {
          const auto& inhits = evt.getValidHandle<StepPointMCCollection>(preSimulatedHitTag_);
          steppingOutput_->reserve(inhits->size());
          for(const auto& hit: *inhits) {
            steppingOutput_->emplace_back(hit);
          }
        }

      }
    }

    void IOHelper::put(art::Event& evt) {
      if(steppingOutput_) {
        evt.put(std::move(steppingOutput_), steppingOutputName_);
      }
    }

    void IOHelper::deleteCutsData(){
      if(steppingOutput_) {
        steppingOutput_ = nullptr;
      }
    }

    void IOHelper::finishConstruction(const CLHEP::Hep3Vector& mu2eOriginInWorld) {
      mu2eOrigin_ = mu2eOriginInWorld;
    }

    void IOHelper::addHit(const G4Step *aStep) {
      if(steppingOutput_->size() < mu2elimits_->maxStepPointCollectionSize()) {

        G4VProcess const* process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
        if(!process) {
          throw cet::exception("GEANT4")<<"ProcessDefinedStep: process not specified for particle "
                                        << aStep->GetTrack()->GetParticleDefinition()->GetParticleName()
                                        << " in file "<<__FILE__<<" line "<<__LINE__
                                        <<" function "<<__func__<<"()\n";
        }

        ProcessCode endCode = ProcessCode::findByName(process->GetProcessName());

        // The point's coordinates are saved in the mu2e coordinate system.
        steppingOutput_->
          push_back(StepPointMC(spHelper_->particlePtr(aStep->GetTrack()),
                                aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(),
                                aStep->GetTotalEnergyDeposit(),
                                aStep->GetNonIonizingEnergyDeposit(),
                                0., // visible energy deposit; used in scintillators
                                aStep->GetPreStepPoint()->GetGlobalTime(),
                                aStep->GetPreStepPoint()->GetProperTime(),
                                aStep->GetPreStepPoint()->GetPosition() - mu2eOrigin_,
                                aStep->GetPostStepPoint()->GetPosition() - mu2eOrigin_,
                                aStep->GetPreStepPoint()->GetMomentum(),
                                aStep->GetPostStepPoint()->GetMomentum(),
                                aStep->GetStepLength(),
                                endCode
                                ));
      }
      else {
        if(!overflowWarningPrinted_) {
          overflowWarningPrinted_ = true;
          mf::LogWarning("G4") << "Maximum number of entries reached in steppingOutput collection "
                               << steppingOutputName_ << ": " << steppingOutput_->size() << endl;
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
      virtual void declareProducts(art::ProducesCollector& pc, art::ConsumesCollector& cc) override;
      virtual void beginEvent(const art::Event& evt, const SimParticleHelper& spHelper) override;
      virtual void put(art::Event&  evt) override;
      virtual void deleteCutsData() override;
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
          if(steppingOutput_) {
            addHit(step);
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

    void Union::declareProducts(art::ProducesCollector& pc, art::ConsumesCollector& cc) {
      IOHelper::declareProducts(pc, cc);
      for(auto& cut: cuts_) {
        cut->declareProducts(pc, cc);
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

    void Union::deleteCutsData(){
      IOHelper::deleteCutsData();
      for(auto& cut: cuts_) {
        cut->deleteCutsData();
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
      virtual void declareProducts(art::ProducesCollector& pc, art::ConsumesCollector& cc) override;
      virtual void beginEvent(const art::Event& evt, const SimParticleHelper& spHelper) override;
      virtual void put(art::Event& evt) override;
      virtual void deleteCutsData() override;
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
      if(result && steppingOutput_) {
        addHit(step);
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

    void Intersection::declareProducts(art::ProducesCollector& pc, art::ConsumesCollector& cc) {
      IOHelper::declareProducts(pc, cc);
      for(auto& cut: cuts_) {
        cut->declareProducts(pc, cc);
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

    void Intersection::deleteCutsData(){
      IOHelper::deleteCutsData();
      for(auto& cut: cuts_) {
        cut->deleteCutsData();
      }
    }


    //================================================================
    class PlaneHelper {
    public:
      explicit PlaneHelper(const fhicl::ParameterSet& pset);
      bool cut_impl(const CLHEP::Hep3Vector& pos);
    private:
      std::array<double,3> normal_;
      double offset_;
    };

    PlaneHelper::PlaneHelper(const fhicl::ParameterSet& pset)
      : offset_()
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

    bool PlaneHelper::cut_impl(const CLHEP::Hep3Vector& pos) {
      const bool result =
        (pos.x()*normal_[0] + pos.y()*normal_[1] + pos.z()*normal_[2] >= offset_);

      return result;
    }

    //================================================================
    class Plane: virtual public IMu2eG4Cut,
                 public IOHelper,
                 private PlaneHelper
    {
    public:
      virtual bool steppingActionCut(const G4Step  *step);
      virtual bool stackingActionCut(const G4Track *trk);

      explicit Plane(const fhicl::ParameterSet& pset, const Mu2eG4ResourceLimits& lim);
    };

    Plane::Plane(const fhicl::ParameterSet& pset, const Mu2eG4ResourceLimits& lim)
      : IOHelper(pset, lim)
      , PlaneHelper(pset)
    {}

    bool Plane::steppingActionCut(const G4Step *step) {
      const CLHEP::Hep3Vector& pos = step->GetPostStepPoint()->GetPosition();
      const bool result = cut_impl(pos);
      if(result && steppingOutput_) {
        addHit(step);
      }
      return result;
    }

    bool Plane::stackingActionCut(const G4Track *trk) {
      const CLHEP::Hep3Vector& pos = trk->GetPosition();
      return cut_impl(pos);
    }

    //================================================================
    // This triggers on steps crossing the plane, not just on the
    // post-step point being behind the plane.  Therefore it is a
    // factor of two slower than Plane, but is sometimes useful
    // for its "virtual detector"-like functionality.

    class ObserverPlane: virtual public IMu2eG4Cut,
                         public IOHelper,
                         private PlaneHelper {
      bool doNotCut_;
    public:
      ObserverPlane(const fhicl::ParameterSet& pset, const Mu2eG4ResourceLimits& lim)
        : IOHelper(pset, lim)
        , PlaneHelper(pset)
        , doNotCut_(pset.get<bool>("doNotCut"))
      {}

      virtual bool stackingActionCut(const G4Track *trk) { return false; }
      virtual bool steppingActionCut(const G4Step  *step);
    };

    bool ObserverPlane::steppingActionCut(const G4Step *step) {
      const bool result =
        cut_impl(step->GetPostStepPoint()->GetPosition())
        && !cut_impl(step->GetPreStepPoint()->GetPosition());

      if(result && steppingOutput_) {
        addHit(step);
      }

      return doNotCut_ ? false : result;
    }

    //================================================================
    class VolumeCut: virtual public IMu2eG4Cut,
                     public IOHelper
    {
    public:
      virtual bool steppingActionCut(const G4Step  *step);
      virtual bool stackingActionCut(const G4Track *trk);

      explicit VolumeCut(const fhicl::ParameterSet& pset, bool negate, const Mu2eG4ResourceLimits& lim);
      virtual void finishConstruction(const CLHEP::Hep3Vector& mu2eOriginInWorld) override;
    private:
      std::vector<std::string> volnames_;
      bool negate_;

      // cached pointers to physical volumes on the list
      typedef std::set<const G4VPhysicalVolume*> KillerVolumesCache;
      KillerVolumesCache killerVolumes_;

      bool cut_impl(const G4Track* trk);
    };

    VolumeCut::VolumeCut(const fhicl::ParameterSet& pset, bool negate, const Mu2eG4ResourceLimits& lim)
      : IOHelper(pset, lim)
      , volnames_(pset.get<std::vector<std::string> >("pars"))
      , negate_(negate)
    {}

    void VolumeCut::finishConstruction(const CLHEP::Hep3Vector& mu2eOriginInWorld) {
      IOHelper::finishConstruction(mu2eOriginInWorld);
      for(const auto& vol: volnames_) {
        killerVolumes_.insert(getPhysicalVolumeOrThrow(vol));
      }
    }

    bool VolumeCut::cut_impl(const G4Track* trk) {
      bool result = false;
      const auto vol = trk->GetVolume();
      // Volume is not defined when we are called from the stacking action.
      // This protection is important for the negated case.
      if(vol) {
        KillerVolumesCache::const_iterator p = killerVolumes_.find(vol);
        result = ( p != killerVolumes_.end() );
        if(negate_) result = !result;
      }
      return result;
    }

    bool VolumeCut::steppingActionCut(const G4Step *step) {
      const bool result = cut_impl(step->GetTrack());
      if(result && steppingOutput_) {
        addHit(step);
      }
      return result;
    }

    bool VolumeCut::stackingActionCut(const G4Track *trk) {
      return cut_impl(trk);
    }

    //================================================================
    class ParticleIdCut: virtual public IMu2eG4Cut,
                         public IOHelper
    {
    public:
      virtual bool steppingActionCut(const G4Step  *step);
      virtual bool stackingActionCut(const G4Track *trk);

      explicit ParticleIdCut(const fhicl::ParameterSet& pset, bool negate, const Mu2eG4ResourceLimits& lim);
    private:
      std::vector<int> pdgIds_;
      bool negate_;
      bool cut_impl(const G4Track* trk);
    };

    ParticleIdCut::ParticleIdCut(const fhicl::ParameterSet& pset, bool negate, const Mu2eG4ResourceLimits& lim)
      : IOHelper(pset, lim)
      , pdgIds_(pset.get<std::vector<int> >("pars"))
      , negate_(negate)
    {
      std::sort(pdgIds_.begin(), pdgIds_.end());
    }

    bool ParticleIdCut::cut_impl(const G4Track* trk) {
      const int id(trk->GetDefinition()->GetPDGEncoding());
      bool found = std::binary_search(pdgIds_.begin(), pdgIds_.end(), id);
      return negate_ ?  !found : found;
    }

    bool ParticleIdCut::steppingActionCut(const G4Step *step) {
      const bool result = cut_impl(step->GetTrack());
      if(result && steppingOutput_) {
        addHit(step);
      }
      return result;
    }

    bool ParticleIdCut::stackingActionCut(const G4Track *trk) {
      return cut_impl(trk);
    }

    //================================================================
    template<class AcceptedCharge>
    class ParticleChargeCut: virtual public IMu2eG4Cut,
                             public IOHelper
    {
    public:
      virtual bool steppingActionCut(const G4Step  *step);
      virtual bool stackingActionCut(const G4Track *trk);

      ParticleChargeCut(const fhicl::ParameterSet& pset, const Mu2eG4ResourceLimits& lim)
        : IOHelper(pset, lim)
      {}

    private:
      GlobalConstantsHandle<ParticleDataTable> pdt_;
      typedef std::map<int,bool> PIDCache;
      PIDCache cache_;
      bool cut_impl(const G4Track* trk);
      AcceptedCharge cut_;
    };

    template<class AcceptedCharge>
    bool ParticleChargeCut<AcceptedCharge>::cut_impl(const G4Track* trk) {
      const int pdgId(trk->GetDefinition()->GetPDGEncoding());
      const auto citer = cache_.find(pdgId);
      if(citer == cache_.end()) {
        ParticleDataTable::maybe_ref info = pdt_->particle(pdgId);
        if(!info.isValid()) {
          throw cet::exception("RUNTIME")<<"ParticleDataTable does onot have information for pdgId = "
                                         << pdgId
                                         << " in file "<<__FILE__<<" line "<<__LINE__
                                         <<" function "<<__func__<<"()\n";
        }
        const double charge = info.ref().charge();
        const bool result = cache_[pdgId] = cut_(charge);
        return result;
      }
      else {
        return citer->second;
      }
    }

    template<class AcceptedCharge>
    bool ParticleChargeCut<AcceptedCharge>::steppingActionCut(const G4Step *step) {
      const bool result = cut_impl(step->GetTrack());
      if(result && steppingOutput_) {
        addHit(step);
      }
      return result;
    }

    template<class AcceptedCharge>
    bool ParticleChargeCut<AcceptedCharge>::stackingActionCut(const G4Track *trk) {
      return cut_impl(trk);
    }

    //----------------------------------------------------------------
    struct AcceptCharged {
      bool operator()(double c) { return std::abs(c) > 0.1; }
    };

    struct AcceptNeutral {
      bool operator()(double c) { return !AcceptCharged()(c); }
    };

    //================================================================
    class KineticEnergy: virtual public IMu2eG4Cut,
                         public IOHelper
    {
    public:
      virtual bool steppingActionCut(const G4Step  *step);
      virtual bool stackingActionCut(const G4Track *trk);

      explicit KineticEnergy(const fhicl::ParameterSet& pset, const Mu2eG4ResourceLimits& lim);
    private:
      double cut_;
      bool cut_impl(const G4Track* trk);
    };

    KineticEnergy::KineticEnergy(const fhicl::ParameterSet& pset, const Mu2eG4ResourceLimits& lim)
      : IOHelper(pset, lim)
      , cut_(pset.get<double>("cut"))
    {}

    bool KineticEnergy::cut_impl(const G4Track* trk) {
      return trk->GetKineticEnergy() < cut_;
    }

    bool KineticEnergy::steppingActionCut(const G4Step *step) {
      const bool result = cut_impl(step->GetTrack());
      if(result && steppingOutput_) {
        addHit(step);
      }
      return result;
    }

    bool KineticEnergy::stackingActionCut(const G4Track *trk) {
      return cut_impl(trk);
    }

    //================================================================
    class GlobalTime: virtual public IMu2eG4Cut,
                      public IOHelper
    {
    public:
      virtual bool steppingActionCut(const G4Step  *step);
      virtual bool stackingActionCut(const G4Track *trk);

      explicit GlobalTime(const fhicl::ParameterSet& pset, const Mu2eG4ResourceLimits& lim);
    private:
      double cut_;
      bool cut_impl(const G4Track* trk);
    };

    GlobalTime::GlobalTime(const fhicl::ParameterSet& pset, const Mu2eG4ResourceLimits& lim)
      : IOHelper(pset, lim)
      , cut_(pset.get<double>("cut"))
    {}

    bool GlobalTime::cut_impl(const G4Track* trk) {
      return trk->GetGlobalTime() > cut_;
    }

    bool GlobalTime::steppingActionCut(const G4Step *step) {
      const bool result = cut_impl(step->GetTrack());
      if(result && steppingOutput_) {
        addHit(step);
      }
      return result;
    }

    bool GlobalTime::stackingActionCut(const G4Track *trk) {
      return cut_impl(trk);
    }

    //================================================================
    class PrimaryOnly: virtual public IMu2eG4Cut,
                       public IOHelper
    {
    public:
      virtual bool steppingActionCut(const G4Step  *step);
      virtual bool stackingActionCut(const G4Track *trk);

      explicit PrimaryOnly(const fhicl::ParameterSet& pset, const Mu2eG4ResourceLimits& lim);
    private:
      bool cut_impl(const G4Track* trk);
    };

    PrimaryOnly::PrimaryOnly(const fhicl::ParameterSet& pset, const Mu2eG4ResourceLimits& lim)
      : IOHelper(pset, lim)
    {}

    bool PrimaryOnly::cut_impl(const G4Track* trk) {
      return (trk->GetParentID() != 0);
    }

    bool PrimaryOnly::steppingActionCut(const G4Step *step) {
      const bool result = cut_impl(step->GetTrack());
      if(result && steppingOutput_) {
        addHit(step);
      }
      return result;
    }

    bool PrimaryOnly::stackingActionCut(const G4Track *trk) {
      return cut_impl(trk);
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
      , value_{pset.get<bool>("value")}
    {}

    Constant::Constant(bool val, const Mu2eG4ResourceLimits& lim) : IOHelper(fhicl::ParameterSet(), lim), value_(val) {}

    bool Constant::steppingActionCut(const G4Step *step) {
      if(steppingOutput_) {
        addHit(step);
      }
      return value_;
    }

    bool Constant::stackingActionCut(const G4Track *) {
      return value_;
    }

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
    if(cuttype == "observerPlane") return make_unique<ObserverPlane>(pset, lim);

    if(cuttype == "inVolume") return make_unique<VolumeCut>(pset, false, lim);
    if(cuttype == "notInVolume") return make_unique<VolumeCut>(pset, true, lim);

    if(cuttype == "pdgId") return make_unique<ParticleIdCut>(pset, false, lim);
    if(cuttype == "notPdgId") return make_unique<ParticleIdCut>(pset, true, lim);

    if(cuttype == "isNeutral") return make_unique<ParticleChargeCut<AcceptNeutral> >(pset, lim);
    if(cuttype == "isCharged") return make_unique<ParticleChargeCut<AcceptCharged> >(pset, lim);

    if(cuttype == "kineticEnergy") return make_unique<KineticEnergy>(pset, lim);
    if(cuttype == "globalTime") return make_unique<GlobalTime>(pset, lim);
    if(cuttype == "primary") return make_unique<PrimaryOnly>(pset, lim);

    if(cuttype == "constant") return make_unique<Constant>(pset, lim);

    throw cet::exception("CONFIG")<< "mu2e::createMu2eG4Cuts(): can not parse pset = "<<pset.to_string()<<"\n";
  }

} // end namespace mu2e
