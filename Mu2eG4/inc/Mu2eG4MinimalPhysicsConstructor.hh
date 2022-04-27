#ifndef Mu2eG4_Mu2eG4MinimalPhysicsConstructor_hh
#define Mu2eG4_Mu2eG4MinimalPhysicsConstructor_hh
//
// A Physics constructor that defines the minimal physics
//
// Original author KLG inspired by Rob Kutschke Minimal physics list
//
#include "Geant4/G4VPhysicsConstructor.hh"

namespace mu2e {

  class  Mu2eG4MinimalPhysicsConstructor: public G4VPhysicsConstructor {

  public:

    Mu2eG4MinimalPhysicsConstructor();

    virtual ~Mu2eG4MinimalPhysicsConstructor() = default;

    void ConstructParticle();

    void ConstructProcess();

  };

} // end namespace mu2e

#endif /* Mu2eG4_Mu2eG4MinimalPhysicsConstructor_hh */
