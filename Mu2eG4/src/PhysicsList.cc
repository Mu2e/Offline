//
// Construct particles; construct and register physics processes.
// This is modeled on:
//   $G4INSTALL/examples/novice/N02/include/ExN02PhysicsList.hh
//    with cvs tag: version 1.12 2008/09/22 16:41:20 maire
//
//
// $Id: PhysicsList.cc,v 1.7 2011/05/18 02:27:18 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:18 $
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>

// Mu2e includes
#include "Mu2eG4/inc/PhysicsList.hh"
#include "Mu2eG4/inc/addStepLimiter.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eReflection/inc/Mu2eReflection.hh"

// G4 includes
#include "globals.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4ionIonisation.hh"

#include "G4Decay.hh"

//#include "G4UserSpecialCuts.hh"

using namespace std;

namespace mu2e{

  PhysicsList::PhysicsList( const SimpleConfig& config):
    G4VUserPhysicsList(),
    _config(&config){

    defaultCutValue = 1.0*cm;
    SetVerboseLevel(1);
  }

  PhysicsList::~PhysicsList(){}

  // Called by the RunManager to define all particles.
  void PhysicsList::ConstructParticle(){
    ConstructBosons();
    ConstructLeptons();
    ConstructMesons();
    ConstructBaryons();
    ConstructAllOthers();
  }

  // Called by the RunManager to define and register processes.
  void PhysicsList::ConstructProcess(){
    AddTransportation();
    ConstructEM();
    ConstructGeneral();

    // Add step limiters to a standard list of particles.
    addStepLimiter();
  }

  // Called by the RunManager to set cuts.
  void PhysicsList::SetCuts(){
    //G4VUserPhysicsList::SetCutsWithDefault method sets
    //the default cut value for all particle types
    SetCutsWithDefault();

    if (verboseLevel>0) DumpCutValuesTable();
  }

  // Methods below here are local methods used to break up the
  // above methods into smaller CLHEP::pieces.

  void PhysicsList::ConstructBosons(){

    // pseudo-particles
    G4Geantino::GeantinoDefinition();
    G4ChargedGeantino::ChargedGeantinoDefinition();

    // gamma
    G4Gamma::GammaDefinition();
  }


  void PhysicsList::ConstructLeptons(){

    // Add tau's if needed.
    G4Electron::ElectronDefinition();
    G4Positron::PositronDefinition();
    G4MuonPlus::MuonPlusDefinition();
    G4MuonMinus::MuonMinusDefinition();
    G4NeutrinoE::NeutrinoEDefinition();
    G4AntiNeutrinoE::AntiNeutrinoEDefinition();
    G4NeutrinoMu::NeutrinoMuDefinition();
    G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
  }


  void PhysicsList::ConstructMesons(){

    // We can add more if neeeded.
    G4PionPlus::PionPlusDefinition();
    G4PionMinus::PionMinusDefinition();
    G4PionZero::PionZeroDefinition();
    G4Eta::EtaDefinition();
    G4EtaPrime::EtaPrimeDefinition();
    G4KaonPlus::KaonPlusDefinition();
    G4KaonMinus::KaonMinusDefinition();
    G4KaonZero::KaonZeroDefinition();
    G4AntiKaonZero::AntiKaonZeroDefinition();
    G4KaonZeroLong::KaonZeroLongDefinition();
    G4KaonZeroShort::KaonZeroShortDefinition();
  }


  void PhysicsList::ConstructBaryons(){

    // We can add more if needed.
    G4Proton::ProtonDefinition();
    G4AntiProton::AntiProtonDefinition();

    G4Neutron::NeutronDefinition();
    G4AntiNeutron::AntiNeutronDefinition();
  }

  void PhysicsList::ConstructAllOthers(){
    // A place holder for nuclei, ions etc.
  }


  // Electromagnetic processes.
  void PhysicsList::ConstructEM(){

    //
    // are we doing Mu2eReflection?

    if (!( _config->get<bool>("mu2eReflection",false))) {
      // Loop over all defined particle types.
      //      cout << "defining physics processes" << endl; assert (2==1);
      theParticleIterator->reset();
      while( (*theParticleIterator)() ){

        // Properties of this particle type.
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();

        // In the following, do the new's leak?

        // Define processes for each particle type.
        if (particleName == "gamma") {

          pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
          pmanager->AddDiscreteProcess(new G4ComptonScattering);
          pmanager->AddDiscreteProcess(new G4GammaConversion);

        } else if (particleName == "e-") {
          pmanager->AddProcess(new G4eMultipleScattering, -1, 1, 1);
          pmanager->AddProcess(new G4eIonisation,         -1, 2, 2);
          pmanager->AddProcess(new G4eBremsstrahlung,     -1, 3, 3);
        } else if (particleName == "e+") {
          pmanager->AddProcess(new G4eMultipleScattering, -1, 1, 1);
          pmanager->AddProcess(new G4eIonisation,         -1, 2, 2);
          pmanager->AddProcess(new G4eBremsstrahlung,     -1, 3, 3);
          pmanager->AddProcess(new G4eplusAnnihilation,    0,-1, 4);

        } else if( particleName == "mu+" ||
                   particleName == "mu-"    ) {
          pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
          pmanager->AddProcess(new G4MuIonisation,        -1, 2, 2);
          pmanager->AddProcess(new G4MuBremsstrahlung,    -1, 3, 3);
          pmanager->AddProcess(new G4MuPairProduction,    -1, 4, 4);

        } else if( particleName == "proton" ||
                   particleName == "pi-" ||
                   particleName == "pi+"    ) {
          pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
          pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);
          pmanager->AddProcess(new G4hBremsstrahlung,     -1, 3, 3);
          pmanager->AddProcess(new G4hPairProduction,     -1, 4, 4);

        } else if( particleName == "alpha" ||
                   particleName == "He3" ||
                   particleName == "GenericIon" ) {
          pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
          pmanager->AddProcess(new G4ionIonisation,       -1, 2, 2);

        } else if ((!particle->IsShortLived()) &&
                   (particle->GetPDGCharge() != 0.0) &&
                   (particle->GetParticleName() != "chargedgeantino")) {
          pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
          pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);
        }
      }
    }

  }

  void PhysicsList::ConstructGeneral(){

    // Does this leak?
    G4Decay* theDecayProcess = new G4Decay();
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      if (theDecayProcess->IsApplicable(*particle)) {
        pmanager ->AddProcess(theDecayProcess);
        // set ordering for PostStepDoIt and AtRestDoIt
        pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
        pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
      }
    }

    if (_config->get<bool>("mu2eReflection",false))
      {
        //mu2e reflection; segregate code for debugging
        Mu2eReflection* theReflectionProcess = new Mu2eReflection(
                                                                  "TargetFoil_",
                                                                  "ToyDSDownstreamVacuum",
                                                                  .001*meter);
        theParticleIterator->reset();
        while( (*theParticleIterator)() ){
          G4ParticleDefinition* particle = theParticleIterator->value();
          G4ProcessManager* pmanager = particle->GetProcessManager();
          if (theReflectionProcess->IsApplicable(*particle)) {
            pmanager ->AddProcess(theReflectionProcess);
            //set ordering for PostStepDoIt
            pmanager ->SetProcessOrdering(theReflectionProcess, idxPostStep);
          }
        }
      }
  }

} // end namespace mu2e
