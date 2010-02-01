//
// Define a minimal physics list. 
// Just transportation; used for debugging geometry.
//
// $Id: MinimalPhysicsList.cc,v 1.1 2010/02/01 00:16:39 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/02/01 00:16:39 $
//
// Original author Rob Kutschke
//

#include "Mu2eG4/inc/MinimalPhysicsList.hh"
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

    // How to do this for one particle type only:
    // Get the definition for a chargedgeantino.
    //G4ParticleDefinition* particle = 
    //G4ParticleTable::GetParticleTable()->FindParticle("chargedgeantino");
    // Add the step limiter process.
    // Must add a userstep limit to any logical volume in which you want this
    // to be active.
    //G4ProcessManager* pmanager = particle->GetProcessManager();
    //pmanager->AddDiscreteProcess(new G4StepLimiter);

    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
      if ( particleName == "e+" ||
	   particleName == "e-" ||
	   particleName == "mu+" ||
	   particleName == "mu-" ||
	   particleName == "chargedgeantino" ){
	pmanager->AddDiscreteProcess(new G4StepLimiter);
      }
    }
  }
  
  void MinimalPhysicsList::SetCuts(){
  }
  
}  // end namespace mu2e
