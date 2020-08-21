#ifndef Mu2eG4_StepLimiterPhysConstructor_hh
#define Mu2eG4_StepLimiterPhysConstructor_hh
//
// A Physics constructor that adds step limiters to some particles.
//
//
// Original author Rob Kutschke
//
#include "G4VPhysicsConstructor.hh"

namespace mu2e {

  class  StepLimiterPhysConstructor: public G4VPhysicsConstructor {

  public:
    StepLimiterPhysConstructor();
    ~StepLimiterPhysConstructor();

    void ConstructParticle();
    void ConstructProcess();

  };

} // end namespace mu2e

#endif /* Mu2eG4_StepLimiterPhysConstructor_hh */
