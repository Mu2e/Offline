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
//
//---------------------------------------------------------------------------
//
// ClassName:   Mu2eG4DecayMuonsWithSpinPhysicsConstructor based on G4DecayPhysics & F05PhysicsList
//              applied to muons/pions/kaons
//
// Author: KLG
//----------------------------------------------------------------------------
//

#include "Mu2eG4/inc/Mu2eG4DecayMuonsWithSpinPhysicsConstructor.hh"

#include "Geant4/G4ParticleDefinition.hh"
#include "Geant4/G4ProcessManager.hh"

#include "Geant4/G4PionDecayMakeSpin.hh"
#include "Geant4/G4DecayWithSpin.hh"

#include "Geant4/G4DecayTable.hh"
#include "Geant4/G4ProcessTable.hh"

#include "Geant4/G4MuonDecayChannelWithSpin.hh"
#include "Geant4/G4MuonRadiativeDecayChannelWithSpin.hh"

#include "Geant4/G4LeptonConstructor.hh"

#include "Geant4/G4MuonPlus.hh"
#include "Geant4/G4MuonMinus.hh"
#include "Geant4/G4PionPlus.hh"
#include "Geant4/G4PionMinus.hh"
#include "Geant4/G4PionZero.hh"
#include "Geant4/G4KaonPlus.hh"
#include "Geant4/G4KaonMinus.hh"
#include "Geant4/G4KaonZeroLong.hh"

// factory
#include "Geant4/G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(Mu2eG4DecayMuonsWithSpinPhysicsConstructor);

Mu2eG4DecayMuonsWithSpinPhysicsConstructor::Mu2eG4DecayMuonsWithSpinPhysicsConstructor(G4int ver)
  :  G4VPhysicsConstructor("Decay"), verbose(ver), wasActivated(false)
{
  fDecayWithSpinProcess = 0;
}

Mu2eG4DecayMuonsWithSpinPhysicsConstructor::Mu2eG4DecayMuonsWithSpinPhysicsConstructor(const G4String& name, G4int ver)
  :  G4VPhysicsConstructor(name), verbose(ver), wasActivated(false)
{
  fDecayWithSpinProcess = 0;
}

Mu2eG4DecayMuonsWithSpinPhysicsConstructor::~Mu2eG4DecayMuonsWithSpinPhysicsConstructor()
{
  delete fDecayWithSpinProcess;
}

void Mu2eG4DecayMuonsWithSpinPhysicsConstructor::ConstructParticle()
{

  if ( 1 < verbose ) {
    G4cout << "Mu2eG4DecayMuonsWithSpinPhysicsConstructor::ConstructParticle invoked" << G4endl;
  }

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();

  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();

  G4DecayTable* MuonPlusDecayTable = new G4DecayTable();
  MuonPlusDecayTable -> Insert(new
                               G4MuonDecayChannelWithSpin("mu+",0.986));
  MuonPlusDecayTable -> Insert(new
                               G4MuonRadiativeDecayChannelWithSpin("mu+",0.014));
  G4MuonPlus::MuonPlusDefinition() -> SetDecayTable(MuonPlusDecayTable);

  G4DecayTable* MuonMinusDecayTable = new G4DecayTable();
  MuonMinusDecayTable -> Insert(new
                                G4MuonDecayChannelWithSpin("mu-",0.986));
  MuonMinusDecayTable -> Insert(new
                                G4MuonRadiativeDecayChannelWithSpin("mu-",0.014));
  G4MuonMinus::MuonMinusDefinition() -> SetDecayTable(MuonMinusDecayTable);
}

void Mu2eG4DecayMuonsWithSpinPhysicsConstructor::ConstructProcess()
{
  if(wasActivated) { return; }
  wasActivated = true;

  // Add Decay With Spin Process
  fDecayWithSpinProcess = new G4DecayWithSpin();

  G4ProcessTable* processTable = G4ProcessTable::GetProcessTable();

  G4VProcess* decay;
  decay = processTable->FindProcess("Decay",G4MuonPlus::MuonPlus());

  G4ProcessManager* fManager;
  fManager = G4MuonPlus::MuonPlus()->GetProcessManager();

  if (fManager) {
    if (decay) fManager->RemoveProcess(decay);
    fManager->AddProcess(fDecayWithSpinProcess);
    if ( 0 < verbose ) {
      G4cout << __func__ << " replacing Decay with DecayWithSpin for mu plus" << G4endl;
    }
    // set ordering for PostStepDoIt and AtRestDoIt
    fManager ->SetProcessOrdering(fDecayWithSpinProcess, idxPostStep);
    fManager ->SetProcessOrdering(fDecayWithSpinProcess, idxAtRest);
  }

  decay = processTable->FindProcess("Decay",G4MuonMinus::MuonMinus());

  fManager = G4MuonMinus::MuonMinus()->GetProcessManager();

  if (fManager) {
    if (decay) fManager->RemoveProcess(decay);
    fManager->AddProcess(fDecayWithSpinProcess);
    if ( 0 < verbose ) {
      G4cout << __func__ << " replacing Decay with DecayWithSpin for mu minus" << G4endl;
    }
    // set ordering for PostStepDoIt and AtRestDoIt
    fManager ->SetProcessOrdering(fDecayWithSpinProcess, idxPostStep);
    fManager ->SetProcessOrdering(fDecayWithSpinProcess, idxAtRest);
  }

  G4PionDecayMakeSpin* poldecay = new G4PionDecayMakeSpin(); // it includes kaons

  // fixme, do this in a loop, extract repeated code

  decay = processTable->FindProcess("Decay",G4PionPlus::PionPlus());

  fManager = G4PionPlus::PionPlus()->GetProcessManager();

  if (fManager) {
    if (decay) fManager->RemoveProcess(decay);
    fManager->AddProcess(poldecay);
    // set ordering for PostStepDoIt and AtRestDoIt
    fManager ->SetProcessOrdering(poldecay, idxPostStep);
    fManager ->SetProcessOrdering(poldecay, idxAtRest);
  }

  decay = processTable->FindProcess("Decay",G4PionMinus::PionMinus());

  fManager = G4PionMinus::PionMinus()->GetProcessManager();

  if (fManager) {
    if (decay) fManager->RemoveProcess(decay);
    fManager->AddProcess(poldecay);
    // set ordering for PostStepDoIt and AtRestDoIt
    fManager ->SetProcessOrdering(poldecay, idxPostStep);
    fManager ->SetProcessOrdering(poldecay, idxAtRest);
  }


  decay = processTable->FindProcess("Decay",G4KaonPlus::KaonPlus());

  fManager = G4KaonPlus::KaonPlus()->GetProcessManager();

  if (fManager) {
    if (decay) fManager->RemoveProcess(decay);
    fManager->AddProcess(poldecay);
    // set ordering for PostStepDoIt and AtRestDoIt
    fManager ->SetProcessOrdering(poldecay, idxPostStep);
    fManager ->SetProcessOrdering(poldecay, idxAtRest);
  }

  decay = processTable->FindProcess("Decay",G4KaonMinus::KaonMinus());

  fManager = G4KaonMinus::KaonMinus()->GetProcessManager();

  if (fManager) {
    if (decay) fManager->RemoveProcess(decay);
    fManager->AddProcess(poldecay);
    // set ordering for PostStepDoIt and AtRestDoIt
    fManager ->SetProcessOrdering(poldecay, idxPostStep);
    fManager ->SetProcessOrdering(poldecay, idxAtRest);
  }

  decay = processTable->FindProcess("Decay",G4KaonZeroLong::KaonZeroLong());

  fManager = G4KaonZeroLong::KaonZeroLong()->GetProcessManager();

  if (fManager) {
    if (decay) fManager->RemoveProcess(decay);
    fManager->AddProcess(poldecay);
    // set ordering for PostStepDoIt and AtRestDoIt
    fManager ->SetProcessOrdering(poldecay, idxPostStep);
    fManager ->SetProcessOrdering(poldecay, idxAtRest);
  }


}
