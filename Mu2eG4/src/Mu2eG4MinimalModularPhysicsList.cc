//
// Define a minimal physics list.
// Just transportation; used for debugging geometry.
//
//
// Original author Rob Kutschke
//

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"

// Geant4 includes
#include "Mu2eG4/inc/Mu2eG4MinimalModularPhysicsList.hh"
#include "Mu2eG4/inc/addStepLimiter.hh"

#include "Geant4/G4ParticleTypes.hh"
#include "Geant4/G4ProcessManager.hh"
#include "Geant4/G4StepLimiter.hh"
#include "Geant4/G4ParticleTable.hh"

namespace mu2e {
  Mu2eG4MinimalModularPhysicsList::Mu2eG4MinimalModularPhysicsList():  G4VModularPhysicsList(){
    defaultCutValue = 1.0*CLHEP::cm;
    SetVerboseLevel(1);
  }

  Mu2eG4MinimalModularPhysicsList::~Mu2eG4MinimalModularPhysicsList(){
  }

  void Mu2eG4MinimalModularPhysicsList::ConstructParticle(){

    G4ChargedGeantino::Definition();
    G4Electron::Definition();
    G4Positron::Definition();
    G4MuonPlus::Definition();
    G4MuonMinus::Definition();
    G4Gamma::Definition();
    G4Proton::Definition();
    G4AntiProton::Definition();
    G4GenericIon::Definition();

  }

  void Mu2eG4MinimalModularPhysicsList::ConstructProcess(){
    AddTransportation();

    // Add step limiters to a standard list of particles.
    addStepLimiter();

  }

  void Mu2eG4MinimalModularPhysicsList::SetCuts(){
  }

}  // end namespace mu2e
