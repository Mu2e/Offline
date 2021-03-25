//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//---------------------------------------------------------------------------
//
// ClassName:   Mu2eG4MinDEDXPhysicsConstructor
//
// Author:      KLG based on G4EmStandardPhysics
//
//----------------------------------------------------------------------------
//
// This class constructs special limited EM physics
//

#include "Mu2eG4/inc/Mu2eG4MinDEDXPhysicsConstructor.hh"
#include "Geant4/G4SystemOfUnits.hh"
#include "Geant4/G4ParticleDefinition.hh"
#include "Geant4/G4EmParameters.hh"
#include "Geant4/G4LossTableManager.hh"

#include "Geant4/G4ComptonScattering.hh"
#include "Geant4/G4GammaConversion.hh"
#include "Geant4/G4PhotoElectricEffect.hh"

#include "Geant4/G4MuBremsstrahlungModel.hh"
#include "Geant4/G4MuPairProductionModel.hh"
#include "Geant4/G4hBremsstrahlungModel.hh"
#include "Geant4/G4hPairProductionModel.hh"

#include "Geant4/G4eIonisation.hh"
#include "Geant4/G4eBremsstrahlung.hh"
#include "Geant4/G4eplusAnnihilation.hh"
#include "Geant4/G4UAtomicDeexcitation.hh"

#include "Geant4/G4MuIonisation.hh"
#include "Geant4/G4MuBremsstrahlung.hh"
#include "Geant4/G4MuPairProduction.hh"
#include "Geant4/G4hBremsstrahlung.hh"
#include "Geant4/G4hPairProduction.hh"

#include "Geant4/G4hIonisation.hh"
#include "Geant4/G4ionIonisation.hh"
#include "Geant4/G4alphaIonisation.hh"

#include "Geant4/G4Gamma.hh"
#include "Geant4/G4Electron.hh"
#include "Geant4/G4Positron.hh"
#include "Geant4/G4MuonPlus.hh"
#include "Geant4/G4MuonMinus.hh"
#include "Geant4/G4PionPlus.hh"
#include "Geant4/G4PionMinus.hh"
#include "Geant4/G4KaonPlus.hh"
#include "Geant4/G4KaonMinus.hh"
#include "Geant4/G4Proton.hh"
#include "Geant4/G4AntiProton.hh"
#include "Geant4/G4Deuteron.hh"
#include "Geant4/G4Triton.hh"
#include "Geant4/G4He3.hh"
#include "Geant4/G4Alpha.hh"
#include "Geant4/G4GenericIon.hh"

#include "Geant4/G4PhysicsListHelper.hh"
#include "Geant4/G4BuilderType.hh"
#include "Geant4/G4EmModelActivator.hh"

// factory
#include "Geant4/G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(Mu2eG4MinDEDXPhysicsConstructor);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Mu2eG4MinDEDXPhysicsConstructor::Mu2eG4MinDEDXPhysicsConstructor(G4int ver)
  : G4VPhysicsConstructor("Mu2eG4MinDEDXPhysicsConstructor"), verbose(ver)
{
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDefaults();
  param->SetVerbose(verbose);
  param->SetLossFluctuations(false); // special case
  SetPhysicsType(bElectromagnetic);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Mu2eG4MinDEDXPhysicsConstructor::~Mu2eG4MinDEDXPhysicsConstructor()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Mu2eG4MinDEDXPhysicsConstructor::ConstructParticle()
{
  // gamma
  G4Gamma::Definition();

  // leptons
  G4Electron::Definition();
  G4Positron::Definition();
  G4MuonPlus::Definition();
  G4MuonMinus::Definition();

  // mesons
  G4PionPlus::Definition();
  G4PionMinus::Definition();
  G4KaonPlus::Definition();
  G4KaonMinus::Definition();

  // barions
  G4Proton::Definition();
  G4AntiProton::Definition();

  // ions
  G4Deuteron::Definition();
  G4Triton::Definition();
  G4He3::Definition();
  G4Alpha::Definition();
  G4GenericIon::Definition();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Mu2eG4MinDEDXPhysicsConstructor::ConstructProcess()
{
  if(verbose > 1) {
    G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
  }
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  // muon & hadron bremsstrahlung and pair production
  G4MuBremsstrahlung* mub = new G4MuBremsstrahlung();
  G4MuPairProduction* mup = new G4MuPairProduction();
  G4hBremsstrahlung* pib = new G4hBremsstrahlung();
  G4hPairProduction* pip = new G4hPairProduction();
  G4hBremsstrahlung* kb = new G4hBremsstrahlung();
  G4hPairProduction* kp = new G4hPairProduction();
  G4hBremsstrahlung* pb = new G4hBremsstrahlung();
  G4hPairProduction* pp = new G4hPairProduction();

  // Add standard EM Processes
  auto myParticleIterator=GetParticleIterator();
  myParticleIterator->reset();
  while( (*myParticleIterator)() ){
    G4ParticleDefinition* particle = myParticleIterator->value();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {

      ph->RegisterProcess(new G4PhotoElectricEffect(), particle);
      ph->RegisterProcess(new G4ComptonScattering(), particle);
      ph->RegisterProcess(new G4GammaConversion(), particle);

    } else if (particleName == "e-") {

      ph->RegisterProcess(new G4eIonisation(), particle);
      ph->RegisterProcess(new G4eBremsstrahlung(), particle);

    } else if (particleName == "e+") {

      ph->RegisterProcess(new G4eIonisation(), particle);
      ph->RegisterProcess(new G4eBremsstrahlung(), particle);
      ph->RegisterProcess(new G4eplusAnnihilation(), particle);

    } else if (particleName == "mu+" ||
               particleName == "mu-"    ) {

      ph->RegisterProcess(new G4MuIonisation(), particle);
      ph->RegisterProcess(mub, particle);
      ph->RegisterProcess(mup, particle);

    } else if (particleName == "alpha" ||
               particleName == "He3") {

      ph->RegisterProcess(new G4ionIonisation(), particle);

    } else if (particleName == "GenericIon") {

      ph->RegisterProcess(new G4ionIonisation(), particle);

    } else if (particleName == "pi+" ||
               particleName == "pi-" ) {

      ph->RegisterProcess(new G4hIonisation(), particle);
      ph->RegisterProcess(pib, particle);
      ph->RegisterProcess(pip, particle);

    } else if (particleName == "kaon+" ||
               particleName == "kaon-" ) {

      ph->RegisterProcess(new G4hIonisation(), particle);
      ph->RegisterProcess(kb, particle);
      ph->RegisterProcess(kp, particle);

    } else if (particleName == "proton" ||
	       particleName == "anti_proton") {

      ph->RegisterProcess(new G4hIonisation(), particle);
      ph->RegisterProcess(pb, particle);
      ph->RegisterProcess(pp, particle);

    } else if (particleName == "B+" ||
	       particleName == "B-" ||
	       particleName == "D+" ||
	       particleName == "D-" ||
	       particleName == "Ds+" ||
	       particleName == "Ds-" ||
               particleName == "anti_He3" ||
               particleName == "anti_alpha" ||
               particleName == "anti_deuteron" ||
               particleName == "anti_lambda_c+" ||
               particleName == "anti_omega-" ||
               particleName == "anti_sigma_c+" ||
               particleName == "anti_sigma_c++" ||
               particleName == "anti_sigma+" ||
               particleName == "anti_sigma-" ||
               particleName == "anti_triton" ||
               particleName == "anti_xi_c+" ||
               particleName == "anti_xi-" ||
	       particleName == "deuteron" ||
	       particleName == "lambda_c+" ||
               particleName == "omega-" ||
               particleName == "sigma_c+" ||
               particleName == "sigma_c++" ||
               particleName == "sigma+" ||
               particleName == "sigma-" ||
               particleName == "tau+" ||
               particleName == "tau-" ||
	       particleName == "triton" ||
               particleName == "xi_c+" ||
               particleName == "xi-" ) {

      ph->RegisterProcess(new G4hIonisation(), particle);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
