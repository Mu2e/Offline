// Pass events with a hit in the specified virtual detector from the same particle
// passing both downstream and upstream, with the specified minimum momentum
//
// Modeled on FilterStepPointReflection by Andrei Gaponenko
// David Brown, 2014

#include <string>
#include <map>
#include <sstream>

// art includes.
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
#include "GeometryService/inc/GeomHandle.hh"

namespace mu2e {

  //================================================================
  class FilterStepPointReflection : public art::EDFilter {
    art::InputTag  input_;
    double cutMomentumMin_;

    // statistics counters
    unsigned numInputEvents_;
    unsigned numPassedEvents_;
    // define reflection virtual detectors.
    std::vector<int> refvids_;
  public:
    explicit FilterStepPointReflection(const fhicl::ParameterSet& pset);
    virtual bool filter(art::Event& event) override;
    virtual void endJob() override;
  };

  //================================================================
  FilterStepPointReflection::FilterStepPointReflection(const fhicl::ParameterSet& pset)
    : art::EDFilter{pset}
    , input_(pset.get<std::string>("input"))
    , cutMomentumMin_(pset.get<double>("cutMomentumMin"))
    , numInputEvents_(0)
    , numPassedEvents_(0)
  {
  // define the reflection plane as the front of the tracker
    refvids_.push_back(VirtualDetectorId::TT_FrontHollow);
    refvids_.push_back(VirtualDetectorId::TT_FrontPA); 
  }

  //================================================================
  bool FilterStepPointReflection::filter(art::Event& event) {
    bool passed(false);
    std::set<art::Ptr<SimParticle>> downstream, upstream;

    auto ih = event.getValidHandle<StepPointMCCollection>(input_);
    for(const auto& hit : *ih) {
    // find hits in the specified virtual detector above the rqeuired momentum
      if(find(refvids_.begin(),refvids_.end(), hit.volumeId()) != refvids_.end()
	 && hit.momentum().mag() > cutMomentumMin_ ){
	art::Ptr<SimParticle> sp = hit.simParticle();
	if(hit.momentum().z() > 0.0){
	  // downstream particle: look for upstream match
	  if(upstream.find(sp) != upstream.end()){
	    passed = true;
	    break;
	  } else // otherwise record this particle
	    downstream.emplace(sp);
	} else {
	  // upstream particle: look for downstream match
	  if(downstream.find(sp) != downstream.end()){
	    passed = true;
	    break;
	  } else // otherwise record this particle
	    upstream.emplace(sp);
	}
      }
    }

    ++numInputEvents_;
    if(passed) { ++numPassedEvents_; }
    return passed;
  }

  //================================================================
  void FilterStepPointReflection::endJob() {
    mf::LogInfo("Summary")
      <<"FilterStepPointReflection_module: passed "
      <<numPassedEvents_<<" / "<<numInputEvents_<<" events\n";
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FilterStepPointReflection);
