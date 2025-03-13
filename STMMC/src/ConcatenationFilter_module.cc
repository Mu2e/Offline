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
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"


namespace mu2e {
  class ConcatenationFilter : public art::EDFilter {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> STMWaveformDigisTag{Name("STMWaveformDigisTag"), Comment("Input tag of StepPointMCs")};
      };
      using Parameters=art::EDFilter::Table<Config>;
      explicit ConcatenationFilter(const Parameters& pset);
      virtual bool filter(art::Event& event) override;
      virtual void endJob() override;
    private:
      art::ProductToken<STMWaveformDigiCollection> STMWaveformDigisToken;
      uint keptEvents = 0, discardedEvents = 0;
  };

  ConcatenationFilter::ConcatenationFilter(const Parameters& conf) :
    art::EDFilter{conf},
    STMWaveformDigisToken(consumes<STMWaveformDigiCollection>(conf().STMWaveformDigisTag())) {};

  bool ConcatenationFilter::filter(art::Event& event) {
    auto const& STMWaveformDigis = event.getProduct(STMWaveformDigisToken);
    // Only keep events that have non-zero size
    if(STMWaveformDigis.size() > 0) {
      keptEvents++;
      return true;
    }
    else {
      discardedEvents++;
      return false;
    };
  };

  void ConcatenationFilter::endJob() {
    mf::LogInfo log("ConcatenationFilter summary");
    log << "=====ConcatenationFilter summary=====\n";
    log << std::left << std::setw(25) << "\tNo. kept events:     " << keptEvents      << "\n";
    log << std::left << std::setw(25) << "\tNo. discarded events:" << discardedEvents << "\n";
    log << "=====================================\n";
  };
}; // namespace mu2e

DEFINE_ART_MODULE(mu2e::ConcatenationFilter)
