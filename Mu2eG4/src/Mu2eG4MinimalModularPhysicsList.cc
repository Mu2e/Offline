//
// Define a minimal physics list.
// Just transportation plus step limiters; used for debugging geometry etc.
//
//
// Original author Rob Kutschke
//

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"

// Geant4 includes
#include "Offline/Mu2eG4/inc/Mu2eG4MinimalModularPhysicsList.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4MinimalPhysicsConstructor.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4StepLimiterPhysicsConstructor.hh"
//#include "Geant4/G4DecayPhysics.hh"

namespace mu2e {
  Mu2eG4MinimalModularPhysicsList::Mu2eG4MinimalModularPhysicsList():
    G4VModularPhysicsList() {
    defaultCutValue = 1.0*CLHEP::cm;
    SetVerboseLevel(1);
    RegisterPhysics( new Mu2eG4MinimalPhysicsConstructor() );
    RegisterPhysics( new Mu2eG4StepLimiterPhysicsConstructor() );
    // An example how to add Decays
    // RegisterPhysics( new G4DecayPhysics(verboseLevel) );
  }

}  // end namespace mu2e
