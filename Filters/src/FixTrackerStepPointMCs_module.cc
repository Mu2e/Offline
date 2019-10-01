////////////////////////////////////////////////////////////////////////
// Class:       FixTrackerStepPointMCs
// Plugin Type: producer (art v2_06_02)
// File:        FixTrackerStepPointMCs_module.cc
//
// This module fixes the duplication of tracker StepPointMCs that occurred in earlier versions of CompressDigiMCs
// before commit 8dbe0d48e
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

#include <memory>

#include "MCDataProducts/inc/StepPointMCCollection.hh"

namespace mu2e {
  class FixTrackerStepPointMCs;
}


class mu2e::FixTrackerStepPointMCs : public art::EDProducer {
public:
  struct Config {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;

    fhicl::Atom<art::InputTag> trackerSteps{Name("trackerSteps"), Comment("Tracker StepPointMC collection to fix")};
  };
  typedef art::EDProducer::Table<Config> Parameters;

  explicit FixTrackerStepPointMCs(const Parameters & conf);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FixTrackerStepPointMCs(FixTrackerStepPointMCs const &) = delete;
  FixTrackerStepPointMCs(FixTrackerStepPointMCs &&) = delete;
  FixTrackerStepPointMCs & operator = (FixTrackerStepPointMCs const &) = delete;
  FixTrackerStepPointMCs & operator = (FixTrackerStepPointMCs &&) = delete;

  // Required functions.
  void produce(art::Event & event) override;

private:

  Config _conf;
};


mu2e::FixTrackerStepPointMCs::FixTrackerStepPointMCs(const Parameters& conf)
  : EDProducer(conf),
    _conf(conf())
{
  // Call appropriate produces<>() functions here.
  produces<StepPointMCCollection>("tracker");
}

void mu2e::FixTrackerStepPointMCs::produce(art::Event & event)
{

  // Implementation of required member function here.
  auto newTrackerSteps = std::unique_ptr<StepPointMCCollection>(new StepPointMCCollection);

  const auto& oldTrackerStepHandle = event.getValidHandle<StepPointMCCollection>(_conf.trackerSteps());
  const auto& oldTrackerSteps = *oldTrackerStepHandle;
  for (const auto& i_oldTrackerStep : oldTrackerSteps) {
    bool already_added = false;
    for (const auto& i_newTrackerStep : *newTrackerSteps) {
      if ( (i_oldTrackerStep.volumeId() == i_newTrackerStep.volumeId()) &&
	   (i_oldTrackerStep.totalEDep() == i_newTrackerStep.totalEDep())
	 ) {
	already_added = true;
      }
    }
    if (!already_added) {
      newTrackerSteps->push_back(i_oldTrackerStep);
    }
  }
  event.put(std::move(newTrackerSteps), "tracker");
}

DEFINE_ART_MODULE(mu2e::FixTrackerStepPointMCs)
