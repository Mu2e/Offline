// Andrei Gaponenko, 2013

#include "Mu2eG4/inc/SimParticlePrimaryHelper.hh"
#include "art/Framework/Principal/Event.h"

namespace mu2e {
  SimParticlePrimaryHelper::SimParticlePrimaryHelper(const art::Event* event,
                                                     const art::ProductID& simProdID,
                                                     const art::EDProductGetter* sim_prod_getter):
    simProdID_(simProdID),
    event_(event),
    simProductGetter_(sim_prod_getter)
  {}

  void SimParticlePrimaryHelper::addEntryFromGenParticle(const art::ValidHandle<GenParticleCollection>& gensHandle, unsigned genId)
  {
    entries_.emplace_back(art::Ptr<GenParticle>(gensHandle, genId),
                          art::Ptr<SimParticle>(simProdID_) );
  }


  void SimParticlePrimaryHelper::addEntryFromSimParticleId(SimParticleCollection::key_type simId)
  {
    art::Handle<GenParticleCollection> gensHandle;
    entries_.emplace_back(art::Ptr<GenParticle>(gensHandle.id()),
                          art::Ptr<SimParticle>(simProdID_,
                                                simId.asUint(),
                                                simProductGetter_));
  }
}
