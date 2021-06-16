// A G4 primary can be created from a GenParticle or from a
// supported object from a previous simulation stage.  This is a
// bookkeeping helper to provide GenParticle and parent Ptrs for
// SimParticles that are primary in the current simulation stage.
// PrimaryGeneratorAction puts information in, and TrackingAction
// uses it when making SimParticles.   The GenParticle Ptrs point
// to an existing collection, while SimParticle Ptrs point to the
// one which is being produced by the current G4 job.
//
// Andrei Gaponenko, 2013, 2021

#ifndef Mu2eG4_inc_SimParticlePrimaryHelper
#define Mu2eG4_inc_SimParticlePrimaryHelper

#include <vector>
#include <variant>

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Handle.h"

#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/SimParticle.hh"

namespace mu2e {

  class StepPointMC;
  class StageParticle;

  class SimParticlePrimaryHelper {
  public:
    // We need to keep the original art::Ptr just for GenParticles.
    // For other cases a Ptr to the new collection will have to be created anyway, so do not bother.
    typedef std::variant<art::Ptr<GenParticle>,
                         const SimParticle*,
                         const StepPointMC*,
                         const StageParticle*>
    InputParticle;


    SimParticlePrimaryHelper(const art::ProductID& simProdID,
                             const art::EDProductGetter* sim_prod_getter);

    art::Ptr<GenParticle> genParticlePtr(int g4TrkID) const;
    art::Ptr<SimParticle> simParticlePrimaryPtr(int g4TrkID) const;


    InputParticle getEntry(int g4TrkID) const;

    template<class T> void addEntry(T t) {
      entries_.emplace_back(t);
    }

  private:

    typedef std::vector<InputParticle> Entries;
    Entries entries_;

    // need these to create art::Ptr to the new SimParticles
    art::ProductID simProdID_;
    const art::EDProductGetter* simProductGetter_;
  };
}

#endif/*Mu2eG4_inc_SimParticlePrimaryHelper*/
