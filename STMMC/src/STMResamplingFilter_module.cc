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
  class STMResamplingFilter : public art::EDFilter {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> StepPointMCsTag{Name("StepPointMCsTag"), Comment("Input tag of StepPointMCs")};
      };
      using Parameters=art::EDFilter::Table<Config>;
      explicit STMResamplingFilter(const Parameters& pset);
      virtual bool filter(art::Event& event) override;
      virtual void endJob() override;
    private:
      art::ProductToken<StepPointMCCollection> StepPointMCsToken;
      uint keptEvents = 0, discardedEvents = 0;
  };

  STMResamplingFilter::STMResamplingFilter(const Parameters& conf) :
    art::EDFilter{conf},
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().StepPointMCsTag())) {};

  bool STMResamplingFilter::filter(art::Event& event) {
    auto const& StepPointMCs = event.getProduct(StepPointMCsToken);
    // Only keep events that have non-zero size
    if(StepPointMCs.size() > 0) {
      keptEvents++;
      return true;
    }
    else {
      discardedEvents++;
      return false;
    };
  };

  void STMResamplingFilter::endJob() {
    mf::LogInfo log("STMResamplingFilter summary");
    log << "No. kept events:      " << keptEvents      << "\n";
    log << "No. discarded events: " << discardedEvents << "\n";
  };
}; // namespace mu2e

DEFINE_ART_MODULE(mu2e::STMResamplingFilter)
