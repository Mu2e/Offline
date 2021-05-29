// Pass events based on the energy depositon and hit locations in CRV
// 
// Yuri Oksuzian, 2019

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
  class FilterCosmicsStage1 : public art::EDFilter {
    typedef std::vector<art::InputTag> InputTags;
    InputTags inputs_;
    double cutEDepMin_;
    double cutEDepMax_;
    double cutZMin_;
    double cutZMax_;

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

      fhicl::Atom<double> cutZMin {
        Name("cutZMin"),
          Comment("The filter passes events if z-pos of the first step in CRV > cutZMin\n"
                  "By default cutZMin=-inf\n"),
          -std::numeric_limits<double>::max()
          };
      
      fhicl::Atom<double> cutZMax {
        Name("cutZMax"),
          Comment("The filter passes events if z-pos of the first step in CRV < cutZMax\n"
                  "By default cutZMax=inf\n"),
          std::numeric_limits<double>::max()
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
    , cutZMin_(conf().cutZMin())
    , cutZMax_(conf().cutZMax())

    , numInputEvents_(0)
    , numPassedEvents_(0)
  {
    for(const auto& i : conf().inputs()) {
      inputs_.emplace_back(i);
    }
  }

  //================================================================
  bool FilterCosmicsStage1::filter(art::Event& event) {
    bool passedE = false;
    bool passedZ = false;
    double totalEDep = 0;
    for(const auto& cn : inputs_) {
      auto ih = event.getValidHandle<StepPointMCCollection>(cn);
      for(const auto& hit : *ih) {
	// skip steps w/o energy deposition
	if( hit.ionizingEdep() < std::numeric_limits<double>::epsilon())
	  continue;
	// Sum up the total energy deposition in CRV
	totalEDep = totalEDep + hit.ionizingEdep();

      }
      std::vector<mu2e::StepPointMC> stepPoints = *ih;
      if (stepPoints.size()){
	// Check if any steps are inside the specified z-region
	if(stepPoints.at(0).position().z() > cutZMin_ && stepPoints.at(0).position().z() < cutZMax_)
	  passedZ = true;
      }
      else
	passedZ = true; // No hits in CRV. Pass z-cut as the default     
    }
    // Select only events with energy deposition in specified range
    if(totalEDep > cutEDepMin_ && totalEDep < cutEDepMax_) 
      passedE = true;

    bool passed=passedE&passedZ;
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

DEFINE_ART_MODULE(mu2e::FilterCosmicsStage1);
