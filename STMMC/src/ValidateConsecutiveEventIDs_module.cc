// Count number of events mixed with Mix.fcl
// Original author: Pawel Plesniak

// stdlib includes
#include <iostream>
#include <string>

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// exception handling
#include "cetlib_except/exception.h"

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
      fhicl::Atom<art::InputTag> stepPointMCsTagEleBeamCat{Name("stepPointMCsTagEleBeamCat"), Comment("Name of EleBeamCat events to search for")};
    };
    typedef art::EDAnalyzer::Table<Config> Parameters;
    explicit CountMixedEvents(const Parameters& conf);
    void analyze(const art::Event& evt) override;
    unsigned int id = 1;
    size_t tmp = 0;

  private:
    Config _conf;
    art::ProductToken<StepPointMCCollection> StepPointMCsTokenEleBeamCat;
  };

  CountMixedEvents::CountMixedEvents(const Parameters& conf)
    : art::EDAnalyzer(conf),
    StepPointMCsTokenEleBeamCat(consumes<StepPointMCCollection>(conf().stepPointMCsTagEleBeamCat())) {};

  void CountMixedEvents::analyze(const art::Event& event) {
    // Get the hits corresponding to the StepPointMCCollection of interest
    auto const& StepPointMCsEleBeamCat = event.getProduct(StepPointMCsTokenEleBeamCat);
    tmp = StepPointMCsEleBeamCat.size();
    if (event.id().event() != (id + 1) && event.id().event() != 1) {
        std::cout << id << ", " << event.id().event() << std::endl;
        throw cet::exception("consecutive", "event ID is not consecutive\n");
    };
    if (event.id().event() != (id + 1))
      std::cout << "Not consecutive event IDs: " << event.id().event() << ", " << id << std::endl;
    id = event.id().event();
    return;
  };
}; // namespace mu2e

DEFINE_ART_MODULE(mu2e::CountMixedEvents)
