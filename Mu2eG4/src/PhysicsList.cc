//
// Define the physics list with G4 for the Mu2e.
//
// $Id: PhysicsList.cc,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//
// Intial implementation is a trivial physics lists for use 
// in debugging the geometry.
//

#include "Mu2eG4/inc/PhysicsList.hh"
#include "G4ParticleTypes.hh"
#include "G4ProcessManager.hh"
#include "G4StepLimiter.hh"
#include "G4ParticleTable.hh"


namespace mu2e {
  PhysicsList::PhysicsList():  G4VUserPhysicsList(){
    defaultCutValue = 1.0*cm;
    SetVerboseLevel(1);
  }

  PhysicsList::~PhysicsList(){
  }

  void PhysicsList::ConstructParticle(){

    G4ChargedGeantino::ChargedGeantinoDefinition();
    G4Electron::ElectronDefinition();
    G4Positron::PositronDefinition();
    G4MuonPlus::MuonPlusDefinition();
    G4MuonMinus::MuonMinusDefinition();
    G4Gamma::GammaDefinition();
    
  }

  void PhysicsList::ConstructProcess(){
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
  
  void PhysicsList::SetCuts(){
  }
  
}  // end namespace mu2e
