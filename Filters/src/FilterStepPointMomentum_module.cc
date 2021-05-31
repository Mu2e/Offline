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

namespace mu2e {

  //================================================================
  class FilterStepPointMomentum : public art::EDFilter {
    typedef std::vector<art::InputTag> InputTags;
    InputTags inputs_;
    double cutMomentumMin_;
    double cutMomentumMax_;

    // statistics counters
    unsigned numInputEvents_;
    unsigned numPassedEvents_;
  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Sequence<art::InputTag> inputs {
        Name("inputs"),
          Comment("Particles and StepPointMCs mentioned in thise collections will be preserved.")
          };

      fhicl::Atom<double> cutMomentumMin {
        Name("cutMomentumMin"),
          Comment("The filter passes events if any of the step points satisties pmag>cutMomentumMin\n"
                  "By default cutMomentumMin=-inf\n"),
          -std::numeric_limits<double>::max()
          };

      fhicl::Atom<double> cutMomentumMax {
        Name("cutMomentumMax"),
          Comment("The filter passes events if any of the step points satisties pmag<cutMomentumMax\n"
                  "By default cutMomentumMax=inf\n"),
          std::numeric_limits<double>::max()
          };

    };

    using Parameters = art::EDFilter::Table<Config>;
    explicit FilterStepPointMomentum(const Parameters& conf);
    virtual bool filter(art::Event& event) override;
    virtual void endJob() override;
  };

  //================================================================
  FilterStepPointMomentum::FilterStepPointMomentum(const Parameters& conf)
    : art::EDFilter{conf}
    , cutMomentumMin_(conf().cutMomentumMin())
    , cutMomentumMax_(conf().cutMomentumMax())
    , numInputEvents_(0)
    , numPassedEvents_(0)
  {
    for(const auto& i : conf().inputs()) {
      inputs_.emplace_back(i);
    }

  }

  //================================================================
  bool FilterStepPointMomentum::filter(art::Event& event) {
    bool passed = false;
    for(const auto& cn : inputs_) {
      auto ih = event.getValidHandle<StepPointMCCollection>(cn);
      for(const auto& hit : *ih) {
        if(hit.momentum().mag() > cutMomentumMin_ && hit.momentum().mag() < cutMomentumMax_) {
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
  void FilterStepPointMomentum::endJob() {
    mf::LogInfo("Summary")
      <<"FilterStepPointMomentum_module: passed "
      <<numPassedEvents_<<" / "<<numInputEvents_<<" events\n";
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FilterStepPointMomentum);
