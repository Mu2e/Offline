// Calculates the efficiency of a simulation stage
//
// Andrew Edmonds, 2020

#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "cetlib_except/exception.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"

#include "MCDataProducts/inc/GenEventCount.hh"
#include "Mu2eUtilities/inc/TriggerResultsNavigator.hh"
#include "MCDataProducts/inc/SimStageEfficiency.hh" // repurposing SimStageEfficiency for the efficiency

namespace mu2e {

  //================================================================
  class SimStageEfficiencyCalculator : public art::EDProducer {

  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<art::InputTag> trigResultsTag{Name("trigResultsTag"), Comment("Input tag for art::TriggerResults")};
      fhicl::Atom<std::string> trigPath{Name("trigPath"), Comment("Specific trigger path we want to calculate the efficiency for")};
      fhicl::Atom<art::InputTag> genCountTag{Name("genCountTag"), Comment("Input tag for mu2e::GenEventCount")};
      fhicl::Atom<bool> printTriggers{Name("printTriggers"), Comment("Switch to turn on the printing of all the trigger names"), false};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit SimStageEfficiencyCalculator(const Parameters& conf);
    virtual void produce(art::Event& event) override;
    virtual void endSubRun(art::SubRun& subrun) override;
    virtual void endRun(art::Run& run) override;

  private:
    Config _conf;
    
    GenEventCount::count_t _passedEvents;
    GenEventCount::count_t _generatedEvents;
  };

  //================================================================
  SimStageEfficiencyCalculator::SimStageEfficiencyCalculator(const Parameters& conf)
    : art::EDProducer{conf},
    _conf(conf()),
    _passedEvents(0),
    _generatedEvents(0)
  {
    produces<mu2e::SimStageEfficiency, art::InRun>();
  }

  //================================================================
  void SimStageEfficiencyCalculator::produce(art::Event& event) {
    
    auto trigResultsH = event.getValidHandle<art::TriggerResults>(_conf.trigResultsTag());
    const art::TriggerResults* trigResults = trigResultsH.product();
    TriggerResultsNavigator tnav(trigResults);
    for (const auto& i_trigPath : tnav.getTrigPaths()) {
      if (_conf.printTriggers()) {
	std::cout << i_trigPath << std::endl;
      }
      if (i_trigPath == _conf.trigPath()) {
	if (tnav.accepted(i_trigPath)) {
	  ++_passedEvents;
	}
	//	std::cout << i_trigPath << ": wasrun = " << tnav.wasrun(i_trigPath) << " accepted = " << tnav.accepted(i_trigPath) << std::endl;
      }
    }
  }

  //================================================================
  void SimStageEfficiencyCalculator::endSubRun(art::SubRun& subrun) {

    auto genCountH = subrun.getValidHandle<mu2e::GenEventCount>(_conf.genCountTag());
    const auto genCount = *genCountH;
    auto genEvents = genCount.count();
    _generatedEvents += genEvents;
  }

  //================================================================
  void SimStageEfficiencyCalculator::endRun(art::Run& run) {

    double efficiency = (double) _passedEvents / _generatedEvents;
    std::cout << "Generated Events = " << _generatedEvents << std::endl;
    std::cout << "Passed Events = " << _passedEvents << std::endl;
    std::cout << "Efficiency = " << efficiency << std::endl;
    
    run.put(std::unique_ptr<SimStageEfficiency>(new SimStageEfficiency(_passedEvents, _generatedEvents)));

    _passedEvents = 0;
    _generatedEvents = 0;
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::SimStageEfficiencyCalculator);
