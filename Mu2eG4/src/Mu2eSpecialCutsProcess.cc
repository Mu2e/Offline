// Original author KL Genser; based on various Geant4 processes
// A Geant4 process to apply special Mu2e cuts;
// for now meant to be a rest discrete process; not along step aspect foreseen for now

#include<limits>

#include "Offline/Mu2eG4/inc/Mu2eSpecialCutsProcess.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4ResourceLimits.hh"

#include "Geant4/G4ios.hh"
#include "Geant4/G4VParticleChange.hh"
#include "Geant4/G4Track.hh"
#include "Geant4/G4Step.hh"

namespace mu2e{

  Mu2eSpecialCutsProcess::Mu2eSpecialCutsProcess(
                                                 const Mu2eG4ResourceLimits& lim,
                                                 const G4String& aName
                                                 )

    : G4VProcess(aName,fUserDefined),
      mu2elimits_(lim)
  {
    SetProcessSubType(0);

    if (verboseLevel>0) {
      G4cout << GetProcessName() << " is created "<< G4endl;
    }
  }

  Mu2eSpecialCutsProcess::~Mu2eSpecialCutsProcess()
  {}

  //  no operation in  AlongStepGPIL for now
  G4double Mu2eSpecialCutsProcess::AlongStepGetPhysicalInteractionLength(
                                                 const G4Track&,
                                                 G4double  ,
                                                 G4double  ,
                                                 G4double& ,
                                                 G4GPILSelection*
                                                 ){ return -1.0; }

  //  no operation in  AlongStepDoIt for now
  G4VParticleChange* Mu2eSpecialCutsProcess::AlongStepDoIt(
                                   const G4Track& ,
                                   const G4Step&
                                   ) { return nullptr; }

  // PostStep (not AtRest)

  G4double Mu2eSpecialCutsProcess::PostStepGetPhysicalInteractionLength(
                                                                        const G4Track& aTrack,
                                                                        G4double, // previousStepSize
                                                                        G4ForceCondition* condition  )
  {
    // condition is set to "Not Forced"; if ommited will allways trigger
    *condition = NotForced;

    G4double proposedStep = std::numeric_limits<double>::max();

    // Maximum number of steps

    if (aTrack.GetCurrentStepNumber() >= static_cast<G4int>(mu2elimits_.maxStepsPerTrack())) {
      proposedStep = 0.;
    }
    if (verboseLevel>0) {
      G4int prec = G4cout.precision(15);
      G4cout << __func__ << " : "
             << GetProcessName()
             << " : current step "
             << aTrack.GetCurrentStepNumber()
             << ", max step limit "
             << static_cast<G4int>(mu2elimits_.maxStepsPerTrack())
             << ", proposed step "
             << std::setw(24) << std::scientific << proposedStep << std::defaultfloat
             << G4endl;
      G4cout.precision(prec);
    }
    return proposedStep;

  }

  G4VParticleChange* Mu2eSpecialCutsProcess::PostStepDoIt( const G4Track& aTrack,
                                                           const G4Step&  )
  //
  // Kill the current particle; deposit its kinetic energy based on the cut used
  //
  {
    if (verboseLevel>0) {
      G4cout << __func__ << " : "
             << GetProcessName()
             << " : current step "
             << aTrack.GetCurrentStepNumber()
             << ", max step limit "
             << static_cast<G4int>(mu2elimits_.maxStepsPerTrack())
             << G4endl;
    }
    aParticleChange.Initialize(aTrack);
    // if the track reached too many steps deposit its energy
    if (aTrack.GetCurrentStepNumber() >=
        static_cast<G4int>(mu2elimits_.maxStepsPerTrack())) {
      aParticleChange.ProposeEnergy(0.) ;
      aParticleChange.ProposeLocalEnergyDeposit(aTrack.GetKineticEnergy()) ;
    }
    // kill the track either way
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    return &aParticleChange;
  }

  // AtRest

  G4double Mu2eSpecialCutsProcess::AtRestGetPhysicalInteractionLength( const G4Track& aTrack,
                                                                       G4ForceCondition* condition
                                                                       )
  {
    // condition is set to "Not Forced"
    *condition = NotForced;

    G4double proposedStep = std::numeric_limits<double>::max();

    // Maximum number of steps

    if (aTrack.GetCurrentStepNumber() > static_cast<G4int>(mu2elimits_.maxStepsPerTrack())) {
      proposedStep = std::numeric_limits<double>::min();
    }
    if (verboseLevel>0) {
      G4int prec = G4cout.precision(15);
      G4cout << __func__ << " : "
             << GetProcessName()
             << " : current step "
             << aTrack.GetCurrentStepNumber()
             << ", max step limit "
             << static_cast<G4int>(mu2elimits_.maxStepsPerTrack())
             << ", proposed step "
             << std::setw(24) << std::scientific << proposedStep << std::defaultfloat
             << G4endl;
      G4cout.precision(prec);
    }
    return proposedStep;
  }

  G4VParticleChange* Mu2eSpecialCutsProcess::AtRestDoIt( const G4Track& aTrack,
                                                         const G4Step&  )
  //
  // Kill the current particle; deposit its kinetic energy based on the cut used
  //
  {
    if (verboseLevel>0) {
      G4cout << __func__ << " : "
             << GetProcessName()
             << " : current step "
             << aTrack.GetCurrentStepNumber()
             << ", max step limit "
             << static_cast<G4int>(mu2elimits_.maxStepsPerTrack())
             << G4endl;
    }
    aParticleChange.Initialize(aTrack);
    // if the track reached too many steps deposit its energy
    if (aTrack.GetCurrentStepNumber() >=
        static_cast<G4int>(mu2elimits_.maxStepsPerTrack())) {
      aParticleChange.ProposeEnergy(0.) ;
      aParticleChange.ProposeLocalEnergyDeposit(aTrack.GetKineticEnergy()) ;
    }
   // kill the track either way
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    return &aParticleChange;
  }
}
