// Special process to record track info before geant4 post step processes occurred
//
// Original author KL Genser
// // C++ includes
// #include <iostream>

// Mu2e includes
#include "Mu2eG4/inc/UserTrackInformation.hh"
#include "Mu2eG4/inc/Mu2eRecorderProcess.hh"

// G4 includes
#include "G4ios.hh"
#include "G4VParticleChange.hh"
#include "G4Track.hh"
#include "G4Step.hh"

namespace mu2e{

  Mu2eRecorderProcess::Mu2eRecorderProcess(const G4String& aName)
    : G4VContinuousProcess(aName,fUserDefined)
  {
    theProcessSubType = 0;
    if (verboseLevel>0) {
      G4cout << GetProcessName() << " is created "<< G4endl;
    }
  }

  Mu2eRecorderProcess::~Mu2eRecorderProcess()
  {}

  G4VParticleChange* Mu2eRecorderProcess::AlongStepDoIt(
                                                        const G4Track& trk,
                                                        const G4Step& step
                                                        )
  {

    pParticleChange->Initialize(trk);
    if (verboseLevel>0) {
      G4int prec = G4cout.precision(15);
      G4cout << __func__
             << " particle "
             << trk.GetParticleDefinition()->GetParticleName()
             << " totE deposit " << std::fixed << trk.GetStep()->GetTotalEnergyDeposit()
             << " NonIonE deposit " << std::fixed << trk.GetStep()->GetNonIonizingEnergyDeposit()
             << " in " << trk.GetVolume()->GetName()
             << " material " << trk.GetMaterial()->GetName()
             << G4endl;
      G4cout << __func__
             << " KE " << step.GetPostStepPoint()->GetKineticEnergy()
             << " momentum direction " << step.GetPostStepPoint()->GetMomentumDirection()
             << G4endl;
      G4cout.precision(prec);
      //      pParticleChange->DumpInfo();
    }

    UserTrackInformation* ti = dynamic_cast<UserTrackInformation*>(trk.GetUserInformation());
    // the two most important lines getting the intermediate info from the G4PostStepPoint
    ti->SetKineticEnergy(step.GetPostStepPoint()->GetKineticEnergy());
    ti->SetMomentumDirection(step.GetPostStepPoint()->GetMomentumDirection());

    return pParticleChange; // return unchanged

  }

  G4double Mu2eRecorderProcess::GetContinuousStepLimit(const G4Track&, // track,
                                                        G4double previousStepSize,
                                                        G4double, // currentMinimumStep
                                                        G4double& // currentSafety
                                                        )
  {
    return DBL_MAX; // do not limit the step
  }

}
