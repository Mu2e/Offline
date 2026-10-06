// Prototype check: print the format version of every CaloDigiCollection in the event (first few events only).
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include <iostream>

namespace mu2e {
  class CaloDigiVersionDump : public art::EDAnalyzer {
    public:
      explicit CaloDigiVersionDump(fhicl::ParameterSet const& pset) : art::EDAnalyzer{pset} {}
      void analyze(art::Event const& event) override {
        if (++n_ > 3) return;
        for (auto const& h : event.getMany<CaloDigiCollection>())
          std::cout << "[CaloDigiVersionDump] event " << event.id().event() << " " << h.provenance()->inputTag()
                    << " size " << h->size() << " format " << (h->empty() ? -1 : h->front().format()) << std::endl;
      }
    private:
      int n_ = 0;
  };
}
DEFINE_ART_MODULE(mu2e::CaloDigiVersionDump)
