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
    virtual ~Mu2eG4StepLimiterPhysicsConstructor() = default;
    Mu2eG4StepLimiterPhysicsConstructor(const Mu2eG4StepLimiterPhysicsConstructor &) = delete;
    Mu2eG4StepLimiterPhysicsConstructor & operator=(const Mu2eG4StepLimiterPhysicsConstructor &) = delete;
    Mu2eG4StepLimiterPhysicsConstructor & operator=( Mu2eG4StepLimiterPhysicsConstructor && ) = delete;

    void ConstructParticle();
    void ConstructProcess();

  };

} // end namespace mu2e

#endif /* Mu2eG4_Mu2eG4StepLimiterPhysicsConstructor_hh */
