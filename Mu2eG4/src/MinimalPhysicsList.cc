//
// Define a minimal physics list.
// Just transportation; used for debugging geometry.
//
// $Id: MinimalPhysicsList.cc,v 1.6 2013/10/25 18:47:09 genser Exp $
// $Author: genser $
// $Date: 2013/10/25 18:47:09 $
//
// Original author Rob Kutschke
//

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"

// Geant4 includes
#include "Mu2eG4/inc/MinimalPhysicsList.hh"
#include "Mu2eG4/inc/addStepLimiter.hh"

#include "G4ParticleTypes.hh"
#include "G4ProcessManager.hh"
#include "G4StepLimiter.hh"
#include "G4ParticleTable.hh"

namespace mu2e {
  MinimalPhysicsList::MinimalPhysicsList():  G4VModularPhysicsList(){
    defaultCutValue = 1.0*CLHEP::cm;
    SetVerboseLevel(1);
  }

  MinimalPhysicsList::~MinimalPhysicsList(){
  }

  void MinimalPhysicsList::ConstructParticle(){

    G4ChargedGeantino::ChargedGeantinoDefinition();
    G4Electron::ElectronDefinition();
    G4Positron::PositronDefinition();
    G4MuonPlus::MuonPlusDefinition();
    G4MuonMinus::MuonMinusDefinition();
    G4Gamma::GammaDefinition();
    G4Proton::Definition();
    G4AntiProton::Definition();
    G4GenericIon::GenericIonDefinition();

  }

  void MinimalPhysicsList::ConstructProcess(){
    AddTransportation();

    // Add step limiters to a standard list of particles.
    addStepLimiter();

  }

  void MinimalPhysicsList::SetCuts(){
  }

}  // end namespace mu2e
