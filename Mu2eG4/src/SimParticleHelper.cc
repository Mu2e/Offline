// Andrei Gaponenko, 2013

#include "Mu2eG4/inc/SimParticleHelper.hh"
#include "Mu2eG4/inc/Mu2eG4Inputs.hh"

#include "art/Framework/Principal/Event.h"
#include "Geant4/G4Track.hh"
#include "Geant4/G4Threading.hh"

namespace mu2e {

  //================================================================
  SimParticleHelper::SimParticleHelper(unsigned simStage,
                                       const Mu2eG4Inputs& inputs,
                                       const art::ProductID& simID,
                                       const art::Event* event,
                                       const art::EDProductGetter* sim_prod_getter)
    : simStage_(simStage)
    , particleNumberOffset_(0u)
    , simID_(simID)
    , event_(event)
    , simProductGetter_(sim_prod_getter)
  {
    // particleNumberOffset_ is set per event, based on the highest number
    // used by SimParticles in previous simulation stages.
    const auto sph = inputs.inputSimParticles(*event);
    if(sph.isValid() && !sph->empty()) {
      particleNumberOffset_ = sph->back().first.asUint();
    }
  }

  //================================================================
  art::Ptr<SimParticle> SimParticleHelper::particlePtr(const G4Track *trk) const {
    return particlePtrFromG4TrackID(trk->GetTrackID());
  }

  art::Ptr<SimParticle> SimParticleHelper::particlePtrFromG4TrackID(int g4TrkID) const {

    return art::Ptr<SimParticle>(simID_,
                                 particleNumberOffset_ + g4TrkID,
                                 simProductGetter_);

  }

  SimParticleCollection::key_type SimParticleHelper::particleKeyFromG4TrackID(int g4TrkID) const {
    return SimParticleCollection::key_type(g4TrkID + particleNumberOffset_);
  }

  const art::EDProductGetter *SimParticleHelper::productGetter() const {
    return simProductGetter_;
  }

  const art::EDProductGetter *SimParticleHelper::otherProductGetter(art::ProductID otherID) const {
    return event_->productGetter(otherID);
  }
}
