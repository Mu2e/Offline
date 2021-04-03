#ifndef Mu2eG4_Mu2eG4StepLimiterPhysicsConstructor_hh
#define Mu2eG4_Mu2eG4StepLimiterPhysicsConstructor_hh
//
// A Physics constructor that adds step limiters to some particles.
//
//
// Original author Rob Kutschke
//
#include "Geant4/G4VPhysicsConstructor.hh"

namespace mu2e {

  class  Mu2eG4StepLimiterPhysicsConstructor: public G4VPhysicsConstructor {

  public:
    Mu2eG4StepLimiterPhysicsConstructor();
    ~Mu2eG4StepLimiterPhysicsConstructor();

    void ConstructParticle();
    void ConstructProcess();

  };

} // end namespace mu2e

#endif /* Mu2eG4_Mu2eG4StepLimiterPhysicsConstructor_hh */
