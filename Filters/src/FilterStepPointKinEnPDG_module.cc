// Pass events with at least one hit satisfying a min momentum cut.
//
// Andrei Gaponenko, 2013

#include <string>
#include <map>
#include <sstream>

// art includes.
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"

namespace mu2e {

  //================================================================
  double getKineticEnergy(const StepPointMC& hit) {
    // unlike generic conditions, MC particle data
    // should not change run-to-run, so static is safe
    // use static for efficiency
    static GlobalConstantsHandle<ParticleDataTable> pdt;

    ParticleDataTable::maybe_ref info = pdt->particle(hit.simParticle()->pdgId());

    if(!info.isValid()) {
      throw cet::exception("MISSINGINFO")<<"No valid PDG info for hit = "<<hit<<"\n";
    }

    const double mass = info.ref().mass();
    return sqrt(hit.momentum().mag2() + std::pow(mass, 2)) - mass;
  }


  //================================================================
  class FilterStepPointKinEnPDG : public art::EDFilter {
    typedef std::vector<art::InputTag> InputTags;
    InputTags inputs_;

    std::vector<double> cutKEmin_;
    std::vector<PDGCode::type> pdgToDrop_;

    // statistics counters
    unsigned numInputEvents_;
    unsigned numPassedEvents_;
  public:
    explicit FilterStepPointKinEnPDG(const fhicl::ParameterSet& pset);
    virtual bool filter(art::Event& event) override;
    virtual void endJob() override;
  };

  //================================================================
  FilterStepPointKinEnPDG::FilterStepPointKinEnPDG(const fhicl::ParameterSet& pset)
    : art::EDFilter{pset}
    , cutKEmin_(pset.get<std::vector<double> >("cutKEmin", std::vector<double>()))
    , numInputEvents_(0)
    , numPassedEvents_(0)
  {
    typedef std::vector<std::string> VS;
    const VS in(pset.get<VS>("inputs"));
    for(const auto& i : in) {
      inputs_.emplace_back(i);
    }
    const auto drop(pset.get<std::vector<int> >("pdgToDrop", std::vector<int>()));
    for(const auto i : drop) {
      pdgToDrop_.emplace_back(PDGCode::type(i));
    }
     if(drop.empty()) {
      throw cet::exception("BADCONFIG")
        <<"FilterStepPointPDG: either pdgToDrop or pdgToKeep must be specified.\n";
    }

  }

  //================================================================
  bool FilterStepPointKinEnPDG::filter(art::Event& event) {
    bool passed = false;
    for(const auto& cn : inputs_) {
      auto ih = event.getValidHandle<StepPointMCCollection>(cn);
      for(const auto& hit : *ih) {
	const PDGCode::type pdgId = hit.simParticle()->pdgId();
	double ke = getKineticEnergy(hit);
	// Pass all the events that contain high momemntum particles that are on the drop PDG list
	for (unsigned i=0; i<pdgToDrop_.size(); i++){
	  if(ke > cutKEmin_[i] && pdgId == pdgToDrop_[i]){
	    passed = true;
	    break;
	  }	 
	}
	// Pass all the events that contain particles that are not on the drop list
	if(std::find(pdgToDrop_.begin(), pdgToDrop_.end(), pdgId) == pdgToDrop_.end()){
	  passed = true;
	  break;
	}
      }
    }

    ++numInputEvents_;
    if(passed) { ++numPassedEvents_; }
    return passed;
  }

  //================================================================
  void FilterStepPointKinEnPDG::endJob() {
    mf::LogInfo("Summary")
      <<"FilterStepPointKinEnPDG_module: passed "
      <<numPassedEvents_<<" / "<<numInputEvents_<<" events\n";
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FilterStepPointKinEnPDG);
