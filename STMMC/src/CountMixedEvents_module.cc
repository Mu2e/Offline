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
      fhicl::Atom<art::InputTag> stepPointMCsTagEleBeamCat{Name("stepPointMCsTagEleBeamCat"), Comment("Name of EleBeamCat events to search for")};
      fhicl::Atom<art::InputTag> stepPointMCsTagMuBeamCat{Name("stepPointMCsTagMuBeamCat"), Comment("Name of MuBeamCat events to search for")};
      fhicl::Atom<art::InputTag> stepPointMCsTagTargetStopsCat{Name("stepPointMCsTagTargetStopsCat"), Comment("Name of TargetStopsCat events to search for")};
    };
    typedef art::EDAnalyzer::Table<Config> Parameters;
    explicit CountMixedEvents(const Parameters& conf);
    void analyze(const art::Event& evt) override;
    void endJob() override;

  private:
    Config _conf;
    art::ProductToken<StepPointMCCollection> StepPointMCsTokenEleBeamCat, StepPointMCsTokenMuBeamCat, StepPointMCsTokenTargetStopsCat;
    int countEleBeamCat = 0, countMuBeamCat = 0, countTargetStopsCat = 0;
  };

  CountMixedEvents::CountMixedEvents(const Parameters& conf)
    : art::EDAnalyzer(conf),
    StepPointMCsTokenEleBeamCat(consumes<StepPointMCCollection>(conf().stepPointMCsTagEleBeamCat())),
    StepPointMCsTokenMuBeamCat(consumes<StepPointMCCollection>(conf().stepPointMCsTagMuBeamCat())),
    StepPointMCsTokenTargetStopsCat(consumes<StepPointMCCollection>(conf().stepPointMCsTagTargetStopsCat())) {};

  void CountMixedEvents::analyze(const art::Event& event) {
    // Get the hits corresponding to the StepPointMCCollection of interest
    auto const& StepPointMCsEleBeamCat = event.getProduct(StepPointMCsTokenEleBeamCat);
    auto const& StepPointMCsMuBeamCat = event.getProduct(StepPointMCsTokenMuBeamCat);
    auto const& StepPointMCsTargetStopsCat = event.getProduct(StepPointMCsTokenTargetStopsCat);
    if (!StepPointMCsEleBeamCat.empty())
      countEleBeamCat++;
    if (!StepPointMCsMuBeamCat.empty())
      countMuBeamCat++;
    if (!StepPointMCsTargetStopsCat.empty())
      countTargetStopsCat++;
    return;
  };

  void CountMixedEvents::endJob() {
    mf::LogInfo log("CountMixedEvents");
    log << "\n====Number of events kept by input dataset====\n";
    log << "EleBeamCat: " << countEleBeamCat << "\n";
    log << "MuBeamCat: " << countMuBeamCat << "\n";
    log << "TargetStopsCat: " << countTargetStopsCat << "\n";
    log << "\n==============================================\n";
    return;
  };
}; // namespace mu2e

DEFINE_ART_MODULE(mu2e::CountMixedEvents)
