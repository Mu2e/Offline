// Special process to record track info before geant4 post step processes occurred
//
// Original author KL Genser
// // C++ includes
// #include <iostream>

// Mu2e includes
#include "Offline/Mu2eG4/inc/Mu2eG4UserTrackInformation.hh"
#include "Offline/Mu2eG4/inc/Mu2eRecorderProcess.hh"

// G4 includes
#include "Geant4/G4ios.hh"
#include "Geant4/G4VParticleChange.hh"
#include "Geant4/G4Track.hh"
#include "Geant4/G4Step.hh"

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
      G4cout << __func__ << " : "
             << GetProcessName()
             << " : particle "
             << trk.GetParticleDefinition()->GetParticleName()
             << " totE deposit " << std::fixed << trk.GetStep()->GetTotalEnergyDeposit()
             << " NonIonE deposit " << std::fixed << trk.GetStep()->GetNonIonizingEnergyDeposit()
             << std::defaultfloat
             << " in " << trk.GetVolume()->GetName()
             << " material " << trk.GetMaterial()->GetName()
             << G4endl;
      G4cout << __func__
             << " Interim KE " << step.GetPostStepPoint()->GetKineticEnergy()
             << " Interim momentum direction " << step.GetPostStepPoint()->GetMomentumDirection()
             << " Interim position " << step.GetPostStepPoint()->GetPosition()
             << " Interim Global Time " << step.GetPostStepPoint()->GetGlobalTime()
             << " Interim Proper Time " << step.GetPostStepPoint()->GetProperTime()
             << G4endl;
      G4cout << __func__
             << " Track Global Time " << trk.GetGlobalTime()
             << " Proper Time " << trk.GetProperTime()
             << G4endl;
      G4cout.precision(prec);
      //      pParticleChange->DumpInfo();
    }

    Mu2eG4UserTrackInformation* ti = dynamic_cast<Mu2eG4UserTrackInformation*>(trk.GetUserInformation());
    // the most important lines getting the intermediate info from the G4PostStepPoint
    ti->SetKineticEnergy(step.GetPostStepPoint()->GetKineticEnergy());
    ti->SetMomentumDirection(step.GetPostStepPoint()->GetMomentumDirection());
    ti->SetGlobalTime(step.GetPostStepPoint()->GetGlobalTime());
    ti->SetProperTime(step.GetPostStepPoint()->GetProperTime());
    ti->SetPosition(step.GetPostStepPoint()->GetPosition());

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
