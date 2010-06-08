//
// Define a minimal physics list. 
// Just transportation; used for debugging geometry.
//
// $Id: MinimalPhysicsList.cc,v 1.3 2010/06/08 21:20:44 genser Exp $
// $Author: genser $
// $Date: 2010/06/08 21:20:44 $
//
// Original author Rob Kutschke
//

#include "Mu2eG4/inc/MinimalPhysicsList.hh"
#include "Mu2eG4/inc/addStepLimiter.hh"

#include "G4ParticleTypes.hh"
#include "G4ProcessManager.hh"
#include "G4StepLimiter.hh"
#include "G4ParticleTable.hh"

namespace mu2e {
  MinimalPhysicsList::MinimalPhysicsList():  G4VUserPhysicsList(){
    defaultCutValue = 1.0*cm;
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

  }

  void MinimalPhysicsList::ConstructProcess(){
    AddTransportation();

    // Add step limiters to a standard list of particles.
    addStepLimiter();

  }
  
  void MinimalPhysicsList::SetCuts(){
  }
  
}  // end namespace mu2e
