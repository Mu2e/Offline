// Andrei Gaponenko, 2013

#include "Mu2eG4/inc/SimParticlePrimaryHelper.hh"
#include "art/Framework/Principal/Event.h"

namespace mu2e {
  SimParticlePrimaryHelper::SimParticlePrimaryHelper(const art::Event* event,
                                                     const art::ProductID& simProdID,
                                                     const art::Handle<GenParticleCollection>& gensHandle,
                                                     const art::EDProductGetter* sim_prod_getter):
    gensHandle_(gensHandle),
    simProdID_(simProdID),
    event_(event),
    simProductGetter_(sim_prod_getter)
  {}
    
    void SimParticlePrimaryHelper::addEntryFromGenParticle(unsigned genId)
    {        
        entries_.emplace_back(art::Ptr<GenParticle>(gensHandle_, genId),
                              art::Ptr<SimParticle>(simProdID_) );
    }
    
    
    void SimParticlePrimaryHelper::addEntryFromStepPointMC(SimParticleCollection::key_type simId)
    {
        entries_.emplace_back(art::Ptr<GenParticle>(gensHandle_.id()),
                              art::Ptr<SimParticle>(simProdID_,
                                                    simId.asUint(),
                                                    simProductGetter_));
    }
}
