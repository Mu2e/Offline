// Count number of events mixed with Mix.fcl
// Original author: Pawel Plesniak

// stdlib includes
#include <iostream>
#include <string>

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// fhicl includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/MCDataProducts/inc/StepPointMC.hh"


namespace mu2e {
  class CountMixedEvents : public art::EDAnalyzer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Atom<art::InputTag> stepPointMCsTag{Name("stepPointMCsTag"), Comment("Name of StepPointMCs to search")};
    };
    typedef art::EDAnalyzer::Table<Config> Parameters;
    explicit CountMixedEvents(const Parameters& conf);
    void analyze(const art::Event& evt) override;
    void endJob() override;

  private:
    Config _conf;
    art::ProductToken<StepPointMCCollection> StepPointMCsToken;
    int count = 0;
  };

  CountMixedEvents::CountMixedEvents(const Parameters& conf)
    : art::EDAnalyzer(conf),
      StepPointMCsToken(consumes<StepPointMCCollection>(conf().stepPointMCsTag())) {};

  void CountMixedEvents::analyze(const art::Event& event) {
    // Get the hits corresponding to the StepPointMCCollection of interest
    auto const& StepPointMCs = event.getProduct(StepPointMCsToken);
    for (auto step : StepPointMCs)
      count++;
    return;
  };

  void CountMixedEvents::endJob() {
    mf::LogInfo log("CountMixedEvents");
    log << "\n==========Data summary==========\n";
    log << "\nNumber of kept events: " << count;
    log << "\n================================\n";
    return;
  };
}; // namespace mu2e

DEFINE_ART_MODULE(mu2e::CountMixedEvents)
