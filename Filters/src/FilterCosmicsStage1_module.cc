// Pass events based on the energy depositon and hit locations in CRV
//
// Yuri Oksuzian, 2019

#include <string>
#include <map>
#include <limits>
#include <sstream>

// art includes.
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/CrvStep.hh"

namespace mu2e {

  //================================================================
  class FilterCosmicsStage1 : public art::EDFilter {
    typedef std::vector<art::InputTag> InputTags;
    InputTags inputs_;
    bool timecut_;
    art::InputTag dsent_;
    double cutEDepMin_;
    double cutEDepMax_;
    double useCrvSteps_;

    // statistics counters
    unsigned numInputEvents_;
    unsigned numPassedEvents_;
  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::OptionalAtom<art::InputTag> DSEnt {
        Name("DSStepPointMCs"),
          Comment("StepPointMCs recorded when the particle enters the DS region.  The time of the earliest DS Step defines the latest time to make the Crv energy sum.  If not specified, no time cut is made when making the energy sum")
          };

      fhicl::Sequence<art::InputTag> inputs {
        Name("inputs"),
          Comment("Particles and StepPointMCs mentioned in thise collections will be preserved.")
          };

      fhicl::Atom<double> cutEDepMin {
        Name("cutEDepMin"),
          Comment("The filter passes events if total energy deposition in CRV > cutEDepMin\n"
                  "By default cutEDepMin=-inf\n"),
          -std::numeric_limits<double>::max()
          };

      fhicl::Atom<double> cutEDepMax {
        Name("cutEDepMax"),
          Comment("The filter passes events if total energy deposition in CRV < cutEDepMax\n"
                  "By default cutEDepMax=inf\n"),
          std::numeric_limits<double>::max()
          };

      fhicl::Atom<bool> useCrvSteps {
        Name("useCrvSteps"),
          Comment("If true then use CRV steps as input, otherwise use StepPointMC\n"),
          false
          };

    };

    using Parameters = art::EDFilter::Table<Config>;
    explicit FilterCosmicsStage1(const Parameters& conf);
    virtual bool filter(art::Event& event) override;
    virtual void endJob() override;
  };

  //================================================================
  FilterCosmicsStage1::FilterCosmicsStage1(const Parameters& conf)
    : art::EDFilter{conf}
    , cutEDepMin_(conf().cutEDepMin())
    , cutEDepMax_(conf().cutEDepMax())
    , useCrvSteps_(conf().useCrvSteps())
    , numInputEvents_(0)
    , numPassedEvents_(0)
  {
    timecut_ = conf().DSEnt(dsent_);

    for(const auto& i : conf().inputs()) {
      inputs_.emplace_back(i);
    }
  }

  //================================================================
  bool FilterCosmicsStage1::filter(art::Event& event) {
    bool passed = false;
    double totalEDep = 0;

    // Optionally only count the Crv hits before the the earliest StepPointMC at the DS entrance
    double dsent_time = std::numeric_limits<float>::max();
    if(timecut_){
      auto dsh = event.getValidHandle<StepPointMCCollection>(dsent_);
      for(const auto& ds : *dsh)
        dsent_time = std::min(dsent_time,ds.time());
    }

    for(const auto& cn : inputs_) {
      if(useCrvSteps_){
        auto ih = event.getValidHandle<CrvStepCollection>(cn);
        for(const auto& hit : *ih)
          if(hit.startTime()<dsent_time)totalEDep += hit.visibleEDep();
      }
      else{
        auto ih = event.getValidHandle<StepPointMCCollection>(cn);
        for(const auto& hit : *ih)
          if(hit.time()<dsent_time)totalEDep += hit.visibleEDep();
      }
    }
    // Select only events with energy deposition in specified range
    if(totalEDep > cutEDepMin_ && totalEDep < cutEDepMax_)
      passed = true;
    ++numInputEvents_;
    if(passed) { ++numPassedEvents_; }
    return passed;
  }

  //================================================================
  void FilterCosmicsStage1::endJob() {
    mf::LogInfo("Summary")
      <<"FilterCosmicsStage1_module: passed "
      <<numPassedEvents_<<" / "<<numInputEvents_<<" events\n";
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FilterCosmicsStage1)
