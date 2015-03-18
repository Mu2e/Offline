// Andrei Gaponenko, 2015

#include "Mu2eG4/inc/Mu2eG4SteppingCuts.hh"

#include <vector>
#include <algorithm>
#include <stdexcept>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h"

//#include "G4Track.hh"
#include "G4Step.hh"

#include "MCDataProducts/inc/ProcessCode.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"
//#include "Mu2eG4/inc/getPhysicalVolumeOrThrow.hh"

namespace mu2e {
  namespace SteppingCuts {

    using namespace std;

    typedef std::vector<fhicl::ParameterSet> PSVector;

    //================================================================
    IOHelper::IOHelper(const std::string& outputName)
      : outputName_(outputName)
      , spHelper_()
    {}

    void IOHelper::declareProducts(art::EDProducer *parent) {
      if(!outputName_.empty()) {
        parent->produces<StepPointMCCollection>(outputName_);
      }
    }

    void IOHelper::beginEvent(const art::Event& evt, const SimParticleHelper& spHelper) {
      spHelper_ = &spHelper;
      if(!outputName_.empty()) {
        output_ = make_unique<StepPointMCCollection>();
        std::cout<<"IOHelper: Creating output collection "<<outputName_
                 <<" ptr = "<<output_.get()
                 <<std::endl;
        // FIXME: handle pre-generated hits
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

      //FIXME: _currentSize += 1;
      //FIXME:
      //FIXME: if ( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      //FIXME:   if( (_currentSize - _sizeLimit)==1 ) {
      //FIXME:     mf::LogWarning("G4") << "Maximum number of particles reached in " 
      //FIXME:                          << SensitiveDetectorName
      //FIXME:                          << ": "
      //FIXME:                          << _currentSize << endl;
      //FIXME:   }
      //FIXME:   return false;
      //FIXME: }

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

    //================================================================
    std::unique_ptr<Union> Union::maybe_instance(const fhicl::ParameterSet& pset) {
      unique_ptr<Union> res;
      if(pset.get<string>("type", "") == "union") {
        res.reset(new Union(pset));
      }
      return res;
    }

    Union::Union(const fhicl::ParameterSet& pset)
      : IOHelper{pset.get<string>("write", "")}
      {
        PSVector pars = pset.get<PSVector>("pars");
        for(const auto& p: pars) {
          cuts_.emplace_back(createCuts(p));
        }
      }


    bool Union::evaluate(const G4Step *step) {
      bool result = false;
      for(const auto& cut : cuts_) {
        if(cut->evaluate(step)) {
          result = true;

          if(output_) {
            addHit(step, ProcessCode::mu2eKillerVolume);
          }

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
    std::unique_ptr<Intersection> Intersection::maybe_instance(const fhicl::ParameterSet& pset) {
      unique_ptr<Intersection> res;
      if(pset.get<string>("type", "") == "intersection") {
        res.reset(new Intersection(pset));
      }
      return res;
    }

    Intersection::Intersection(const fhicl::ParameterSet& pset)
      : IOHelper{pset.get<string>("write", "")}
      {
        PSVector pars = pset.get<PSVector>("pars");
        for(const auto& p: pars) {
          cuts_.emplace_back(createCuts(p));
        }
      }

    bool Intersection::evaluate(const G4Step *step) {
      bool result = false;
      for(const auto& cut : cuts_) {
        if(cut->evaluate(step)) {
          result = true;

          if(output_) {
            addHit(step, ProcessCode::mu2eKillerVolume);
          }

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
    std::unique_ptr<Plane> Plane::maybe_instance(const fhicl::ParameterSet& pset) {
      unique_ptr<Plane> res;
      if(pset.get<string>("type", "") == "plane") {
        res.reset(new Plane(pset));
      }
      return res;
    }

    Plane::Plane(const fhicl::ParameterSet& pset)
      : IOHelper{pset.get<string>("write", "")}
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
        throw std::runtime_error("SteppingCut::Plane(): point should be a vector of 3 doubles. Error in pset = "+pset.to_string());
      }
      // the cut:
      //             (x-x0)*normal >= 0
      // rewrite as
      //
      //              x*normal >= x0*normal =: offset

      offset_ = std::inner_product(normal_.begin(), normal_.end(), x0.begin(), 0.);
    }

    bool Plane::evaluate(const G4Step *step) {
      const CLHEP::Hep3Vector& pos2 = step->GetPostStepPoint()->GetPosition();

      const bool result =
        (pos2.x()*normal_[0] + pos2.y()*normal_[1] + pos2.z()*normal_[2] >= offset_);

      //std::cout<<"Plane["<<outputName_<<"]::evaluate() = "<<result<<" for pos2 = "<<pos2<<std::endl;

      if(result && output_) {
        addHit(step, ProcessCode::mu2eKillerVolume);
      }

      return result;
    }

    //================================================================
    std::unique_ptr<Constant> Constant::maybe_instance(const fhicl::ParameterSet& pset) {
      unique_ptr<Constant> res;
      if(pset.get<string>("type", "") == "constant") {
        //res = make_unique<Constant>(pset);
        res.reset(new Constant(pset));
      }
      return res;
    }

    Constant::Constant(const fhicl::ParameterSet& pset)
      : IOHelper{pset.get<string>("write", "")}
      , value_{pset.get<double>("value")}
    {}

    Constant::Constant(bool val) : IOHelper{std::string()}, value_(val) {}

    bool Constant::evaluate(const G4Step *step) {
      if(output_) {
        addHit(step, ProcessCode::mu2eKillerVolume);
      }
      return value_;
    }

    //================================================================
    std::unique_ptr<IMu2eG4SteppingCut> createCuts(const fhicl::ParameterSet& pset) {
      std::unique_ptr<IMu2eG4SteppingCut> res;

      if(pset.is_empty()) {
        res = make_unique<Constant>(false); // no cuts
        return res;
      }

      res = Union::maybe_instance(pset);
      if(res) return res;

      res = Intersection::maybe_instance(pset);
      if(res) return res;

      res = Plane::maybe_instance(pset);
      if(res) return res;

      res = Constant::maybe_instance(pset);
      if(res) return res;

      throw std::runtime_error("mu2e::SteppingCuts::createCuts(pset): can not parse pset = "+pset.to_string());
    }


    //================================================================

  }
}
