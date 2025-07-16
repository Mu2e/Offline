#include "Offline/Mu2eG4/inc/Mu2eG4BiasedRDPhysics.hh"

#include "Geant4/G4Radioactivation.hh"
#include "Geant4/G4GenericIon.hh"
#include "Geant4/globals.hh"
#include "Geant4/G4PhysicsListHelper.hh"
#include "Geant4/G4EmParameters.hh"
#include "Geant4/G4VAtomDeexcitation.hh"
#include "Geant4/G4UAtomicDeexcitation.hh"
#include "Geant4/G4LossTableManager.hh"
#include "Geant4/G4NuclearLevelData.hh"
#include "Geant4/G4DeexPrecoParameters.hh"
#include "Geant4/G4NuclideTable.hh"
#include "Geant4/G4PhysicsConstructorFactory.hh"
#include "Geant4/G4SystemOfUnits.hh"

namespace mu2e{

G4_DECLARE_PHYSCONSTR_FACTORY(Mu2eG4BiasedRDPhysics);


Mu2eG4BiasedRDPhysics::Mu2eG4BiasedRDPhysics(const Mu2eG4Config::Physics* phys,
                                             G4int verbose)
 : G4VPhysicsConstructor("G4Radioactivation"),
   phys_{phys},
   verbose_{verbose}
{
  G4EmParameters::Instance()->AddPhysics("World","G4Radioactivation");
  G4DeexPrecoParameters* deex = G4NuclearLevelData::GetInstance()->GetParameters();
  deex->SetStoreICLevelData(true);
  deex->SetStoreAllLevels(true);
  deex->SetInternalConversionFlag(true);
  deex->SetIsomerProduction(true);
  deex->SetCorrelatedGamma(false);
  deex->SetMaxLifeTime(G4NuclideTable::GetInstance()->GetThresholdOfHalfLife()/std::log(2.));

  const G4double meanLife = 1*picosecond;
  G4NuclideTable::GetInstance()->SetMeanLifeThreshold(meanLife);
  G4NuclideTable::GetInstance()->SetLevelTolerance(1.0*eV);

  // define flags for the atomic de-excitation module
  G4EmParameters::Instance()->SetAugerCascade(true);
  G4EmParameters::Instance()->SetDeexcitationIgnoreCut(true);
}


Mu2eG4BiasedRDPhysics::Mu2eG4BiasedRDPhysics()
 : Mu2eG4BiasedRDPhysics(nullptr)
{ }


void Mu2eG4BiasedRDPhysics::ConstructParticle()
{
  G4GenericIon::GenericIon();
}

void Mu2eG4BiasedRDPhysics::ConstructProcess()
{
  G4LossTableManager* man = G4LossTableManager::Instance();
  G4VAtomDeexcitation* ad = man->AtomDeexcitation();
  if (!ad) {
    G4EmParameters::Instance()->SetAugerCascade(true);
    ad = new G4UAtomicDeexcitation();
    man->SetAtomDeexcitation(ad);
    ad->InitialiseAtomicDeexcitation();
  }

  //ownership ofobject is trasnferred to Geant4
  auto process = new G4Radioactivation();
  process->SetVerboseLevel(verbose_);
  process->SetARM(false); //Atomic Rearangement

   if (phys_){
     process->SetHLThreshold(phys_->radiationHLT());
     process->SetBRBias(phys_->radiationBRBias());
     process->SetSplitNuclei(phys_->radiationNsplit());
     process->SetSourceTimeProfile(phys_->beamTimeProfileRad());
     process->SetDecayBias(phys_->coolTimeProfileRad());
  }

  G4PhysicsListHelper::GetPhysicsListHelper()->
    RegisterProcess(process, G4GenericIon::GenericIon());
}


}
