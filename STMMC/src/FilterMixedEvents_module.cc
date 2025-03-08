// Removes empty StepPointMC collections with the associated tag
// Original author: Pawel Plesniak

// stdlib includes
#include <iostream>
#include <string>

// art includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// fhicl includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/MCDataProducts/inc/StepPointMC.hh"


namespace mu2e {
    class FilterMixedEvents : public art::EDFilter {
    public:
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        struct Config {
            fhicl::Atom<art::InputTag> stepPointMCsTagEleBeamCat{Name("stepPointMCsTagEleBeamCat"), Comment("Name of EleBeamCat events to search for")};
            fhicl::Atom<art::InputTag> stepPointMCsTagMuBeamCat{Name("stepPointMCsTagMuBeamCat"), Comment("Name of MuBeamCat events to search for")};
            fhicl::Atom<art::InputTag> stepPointMCsTagTargetStopsCat{Name("stepPointMCsTagTargetStopsCat"), Comment("Name of TargetStopsCat events to search for")};
        };
        using Parameters=art::EDFilter::Table<Config>;
        explicit FilterMixedEvents(const Parameters& pset);
        virtual bool filter(art::Event& event) override;
        virtual void endJob() override;
    private:
        art::ProductToken<StepPointMCCollection> StepPointMCsTokenEleBeamCat, StepPointMCsTokenMuBeamCat, StepPointMCsTokenTargetStopsCat;
        uint keptEvents = 0, discardedEvents = 0;
  };

  FilterMixedEvents::FilterMixedEvents(const Parameters& conf) :
    art::EDFilter{conf},
    StepPointMCsTokenEleBeamCat(consumes<StepPointMCCollection>(conf().stepPointMCsTagEleBeamCat())),
    StepPointMCsTokenMuBeamCat(consumes<StepPointMCCollection>(conf().stepPointMCsTagMuBeamCat())),
    StepPointMCsTokenTargetStopsCat(consumes<StepPointMCCollection>(conf().stepPointMCsTagTargetStopsCat())) {};

  bool FilterMixedEvents::filter(art::Event& event) {
    // Get the hits corresponding to the StepPointMCCollection of interest
    auto const& StepPointMCsEleBeamCat = event.getProduct(StepPointMCsTokenEleBeamCat);
    auto const& StepPointMCsMuBeamCat = event.getProduct(StepPointMCsTokenMuBeamCat);
    auto const& StepPointMCsTargetStopsCat = event.getProduct(StepPointMCsTokenTargetStopsCat);
    // Only keep events that have non-zero size
    if(StepPointMCsEleBeamCat.empty() && StepPointMCsMuBeamCat.empty() && StepPointMCsTargetStopsCat.empty()) {
        discardedEvents++;
        return false;
    }
    else {
        keptEvents++;
        return true;
    };
  };

  void FilterMixedEvents::endJob() {
    mf::LogInfo log("FilterMixedEvents summary");
    log << "==========================================\n";
    log << "No. kept events:      " << keptEvents      << "\n";
    log << "No. discarded events: " << discardedEvents << "\n";
    log << "==========================================\n";
  };
}; // namespace mu2e

DEFINE_ART_MODULE(mu2e::FilterMixedEvents)
