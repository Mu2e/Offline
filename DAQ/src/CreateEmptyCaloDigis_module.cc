// ======================================================================
//
//  Dummy module to create an empty CaloDigi collection when the raw data doesn't have one.
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include <memory>

// ======================================================================
namespace art {
  class CreateEmptyCaloDigis;
}

class art::CreateEmptyCaloDigis : public EDProducer {
  public:
    struct Config {
    };
    explicit CreateEmptyCaloDigis(const art::EDProducer::Table<Config>& config);
    virtual ~CreateEmptyCaloDigis() {}
    void produce(Event&) override;
};
art::CreateEmptyCaloDigis::CreateEmptyCaloDigis(const art::EDProducer::Table<Config>& config) :
 EDProducer{config} {
  produces<mu2e::CaloDigiCollection>();
 }
// ----------------------------------------------------------------------
void art::CreateEmptyCaloDigis::produce(Event& event) {
  std::unique_ptr<mu2e::CaloDigiCollection> calo_digis(new mu2e::CaloDigiCollection);
  event.put(std::move(calo_digis));
}
DEFINE_ART_MODULE(art::CreateEmptyCaloDigis)
