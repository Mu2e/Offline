// A G4 primary can be created from a GenParticle or from a
// StepPointMC from a previous simulation stage.  This is a
// bookkeeping helper to provide GenParticle and parent Ptrs for
// SimParticles that are primary in the current simulation stage.
// PrimaryGeneratorAction puts information in, and TrackingAction
// uses it when making SimParticles.   The GenParticle Ptrs point
// to an existing collection, while SimParticle Ptrs point to the
// one which is being produced by the current G4 job.
//
// Andrei Gaponenko, 2013

#ifndef Mu2eG4_inc_SimParticlePrimaryHelper
#define Mu2eG4_inc_SimParticlePrimaryHelper

#include <vector>

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Handle.h"

#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

namespace art { class Event; }

namespace mu2e {

  class SimParticlePrimaryHelper {
  public:

    SimParticlePrimaryHelper(const art::Event* event,
                             const art::ProductID& simProdID,
                             const art::EDProductGetter* sim_prod_getter);

    unsigned numPrimaries() const { return entries_.size(); }

    const art::Ptr<GenParticle>& genParticlePtr(int g4TrkID) const {
      return entries_.at(g4TrkID - 1).genParticlePtr;
    }

    const art::Ptr<SimParticle>& simParticlePrimaryPtr(int g4TrkID) const {
      return entries_.at(g4TrkID - 1).simParticlePrimaryPtr;
    }

    void addEntryFromGenParticle(const art::ValidHandle<GenParticleCollection>& gensHandle, unsigned genId);

    void addEntryFromSimParticleId (SimParticleCollection::key_type simId);

  private:

    struct Entry {
      art::Ptr<GenParticle> genParticlePtr;
      art::Ptr<SimParticle> simParticlePrimaryPtr;
      Entry(const art::Ptr<GenParticle>& g, const art::Ptr<SimParticle>& p)
        : genParticlePtr(g), simParticlePrimaryPtr(p)
      {}
    };

    typedef std::vector<Entry> Entries;
    Entries entries_;

    // need these to create art::Ptr to SimParticles
    art::ProductID simProdID_;
    const art::Event* event_;
    const art::EDProductGetter* simProductGetter_;

  };
}

#endif/*Mu2eG4_inc_SimParticlePrimaryHelper*/
