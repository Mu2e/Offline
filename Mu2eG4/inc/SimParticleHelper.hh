// The job of this class is to help with SimParticle creation.
// It translates between G4 track IDs and keys of SimParticles
// produced by the current job.  It also keeps the simulation
// stage number.
//
// Andrei Gaponenko, 2013

#ifndef Mu2eG4_inc_SimParticleHelper
#define Mu2eG4_inc_SimParticleHelper

#include "canvas/Persistency/Common/Ptr.h"

#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

class G4Track;
namespace art { class Event; }
namespace art { class ProductID; }
namespace mu2e { class Mu2eG4Inputs; }

namespace mu2e {
  class SimParticleHelper {
    unsigned simStage_;
    unsigned particleNumberOffset_;
    art::ProductID simID_;
    const art::Event* event_;
    const art::EDProductGetter* simProductGetter_;

  public:
    SimParticleHelper(unsigned simStage,
                      const Mu2eG4Inputs& inputs,
                      const art::ProductID& simID,
                      const art::Event* event,
                      const art::EDProductGetter* sim_prod_getter);

    art::Ptr<SimParticle> particlePtr(const G4Track* trk) const;

    art::Ptr<SimParticle> particlePtrFromG4TrackID(int g4TrkID) const;

    SimParticleCollection::key_type particleKeyFromG4TrackID(int g4TrkID) const;

    // of the SimParticleCollection the current G4 run produces
    art::ProductID productID() const { return simID_; }
    // of the SimParticleCollection the current G4 run produces
    const art::EDProductGetter *productGetter() const;

    // For arbitrary product in this event - e.g. the "old" SimParticleCollection
    const art::EDProductGetter *otherProductGetter(art::ProductID otherID) const;

    unsigned simStage() const { return simStage_; }
  };
}

#endif/*Mu2eG4_inc_SimParticleHelper*/
