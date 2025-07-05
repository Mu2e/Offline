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

// Mu2e includes
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
// Framework includes
#include "canvas/Persistency/Common/Ptr.h"

namespace art { template <typename T> class Ptr; }

namespace mu2e {

  class SimParticle;

  class SimParticleGetTau {
  public:
    typedef std::vector<StepPointMCCollection> VspMC;

    static double calculate( const art::Ptr<SimParticle>& p,
                             const VspMC& hitColls,
                             const std::vector<int>& decayOffCodes = std::vector<int>(),
                             const PhysicsParams& gc = *GlobalConstantsHandle<PhysicsParams>());

    static double calculate( const StepPointMC& sp,
                             const VspMC& hitColls,
                             const std::vector<int>& decayOffCodes = std::vector<int>(),
                             const PhysicsParams& gc = *GlobalConstantsHandle<PhysicsParams>());

    // using a sim particle collection from a multi-stage job
    static double calculate(const art::Ptr<SimParticle>& p,
                            const std::vector<int>& decayOffCodes = std::vector<int>(),
                            const PhysicsParams& gc = *GlobalConstantsHandle<PhysicsParams>());
  private:
    static double getMultiStageTau( const art::Ptr<SimParticle>& sp,
                                    const VspMC& hitColls,
                                    const std::vector<int>& decayOffCodes,
                                    const PhysicsParams& gc );
  };
}

#endif/*Mu2eUtilities_SimParticleGetTau_hh*/
