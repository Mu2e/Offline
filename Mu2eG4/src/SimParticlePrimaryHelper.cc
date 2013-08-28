// Andrei Gaponenko, 2013

#include "Mu2eG4/inc/SimParticlePrimaryHelper.hh"
#include "art/Framework/Principal/Event.h"

namespace mu2e {
  SimParticlePrimaryHelper::SimParticlePrimaryHelper(const art::Event& event,
                                                     const art::ProductID& simProdID,
                                                     const art::Handle<GenParticleCollection>& gensHandle)
    : gensHandle_(gensHandle)
    , simProdID_(simProdID)
    , event_(&event)
  {}


  void SimParticlePrimaryHelper::addEntry(unsigned genId,
                                          SimParticleCollection::key_type simId)
  {
    entries_.emplace_back(
                          // Note that art::Ptr is supposed to have a valid product ID to be a NULL Ptr.
                          (genId == -1u) ? art::Ptr<GenParticle>(gensHandle_.id())
                          : art::Ptr<GenParticle>(gensHandle_, genId),

                          (simId == SimParticleCollection::key_type()) ? art::Ptr<SimParticle>(simProdID_)
                          : art::Ptr<SimParticle>(simProdID_, simId.asUint(), event_->productGetter(simProdID_))
                          );
  }
}
