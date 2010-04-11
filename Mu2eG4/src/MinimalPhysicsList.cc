//
// Define a minimal physics list. 
// Just transportation; used for debugging geometry.
//
// $Id: MinimalPhysicsList.cc,v 1.2 2010/04/11 15:16:55 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/11 15:16:55 $
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

  }

  void MinimalPhysicsList::ConstructProcess(){
    AddTransportation();

    // Add step limiters to a standard list of particles.
    addStepLimiter();

  }
  
  void MinimalPhysicsList::SetCuts(){
  }
  
}  // end namespace mu2e
