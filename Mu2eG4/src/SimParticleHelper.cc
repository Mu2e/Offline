// Andrei Gaponenko, 2013

#include "Mu2eG4/inc/SimParticleHelper.hh"

#include "art/Framework/Principal/Event.h"
#include "G4Track.hh"

namespace mu2e {
  SimParticleHelper::SimParticleHelper(unsigned particleNumberOffset,
                                       const art::ProductID& simID,
                                       const art::Event& event)
    : particleNumberOffset_(particleNumberOffset)
    , simID_(simID)
    , event_(&event)
  {}

  art::Ptr<SimParticle> SimParticleHelper::particlePtr(const G4Track *trk) const {
    return particlePtrFromG4TrackID(trk->GetTrackID());
  }

  art::Ptr<SimParticle> SimParticleHelper::particlePtrFromG4TrackID(int g4TrkID) const {
    return art::Ptr<SimParticle>(simID_,
                                 particleNumberOffset_ + g4TrkID,
                                 event_->productGetter(simID_));
  }

  SimParticleCollection::key_type SimParticleHelper::particleKeyFromG4TrackID(int g4TrkID) const {
    return SimParticleCollection::key_type(g4TrkID + particleNumberOffset_);
  }

  const art::EDProductGetter *SimParticleHelper::productGetter() const {
    return event_->productGetter(simID_);
  }

  const art::EDProductGetter *SimParticleHelper::otherProductGetter(art::ProductID otherID) const {
    return event_->productGetter(otherID);
  }
}
