// static member function to determine proper lifetime of particles
// from multi-staged simulations.  Requires the SimParticle instance
// as the first argument, and then a set of StepPointMC collections
// corresponding to the various stages of the simulation.
//
// Original authors: Kyle Knoepfel
//                   Andrei Gaponenko

#ifndef Mu2eUtilities_SimParticleGetTau_hh
#define Mu2eUtilities_SimParticleGetTau_hh

// C++ includes
#include <vector>

// Framework includes
#include "art/Persistency/Common/Ptr.h"

// Mu2e includes
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/PhysicsParams.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

namespace mu2e {

  class SimParticle;

  class SimParticleGetTau {
    typedef std::vector<StepPointMCCollection> VspMC;
    
  public:
    static double calculate( const art::Ptr<SimParticle>& p, 
                             const VspMC& hitColls, 
                             const std::vector<int>& decayOffCodes = std::vector<int>(),
                             const PhysicsParams& gc = *GlobalConstantsHandle<PhysicsParams>());
  };
}

#endif/*Mu2eUtilities_SimParticleGetTau_hh*/
