//
//  Trivial module to insert an empty PrimaryParticle object in the event, used for NoPrimary production
//
// mu2e
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
namespace mu2e {
  class NullMCPrimary : public art::EDProducer {
    public:
      struct Config {
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit NullMCPrimary(const Parameters& conf);
      void produce(art::Event& evt) override;
    private:
  };

  NullMCPrimary::NullMCPrimary(const Parameters& config )  :
    art::EDProducer{config}
    {
      produces <PrimaryParticle>();
    }

  void NullMCPrimary::produce(art::Event& event) {
    // create empty output object
    PrimaryParticle pp;
    event.put(std::make_unique<PrimaryParticle>(pp));
  }
}
DEFINE_ART_MODULE(mu2e::NullMCPrimary)
