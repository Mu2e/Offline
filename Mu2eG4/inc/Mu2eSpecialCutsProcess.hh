#ifndef Mu2eG4_Mu2eSpecialCutsProcess_hh
#define Mu2eG4_Mu2eSpecialCutsProcess_hh

#include "Offline/Mu2eG4/inc/Mu2eG4ResourceLimits.hh"

#include "Geant4/G4VProcess.hh"


// A process to apply various step related cuts and optionally deposit
// particle kinetic energy; added in toggleProcesses addUserProcesses
// invoked from Mu2eG4CustomizationPhysicsConstructor


// Original author KL Genser; based on various Geant4 processes

namespace mu2e {

  class Mu2eSpecialCutsProcess : public G4VProcess

  {

  public:

    explicit Mu2eSpecialCutsProcess(
                                    const Mu2eG4ResourceLimits& lim,
                                    const G4String& processName ="mu2eSpecialCutsProcess"
                                    );

    virtual ~Mu2eSpecialCutsProcess();

    // virtual G4bool IsApplicable(const G4ParticleDefinition& aParticleType) override;
    // not implemented, i.e., always true however see that is is not
    // added to stable AtRest process vectors

    virtual inline G4double PostStepGetPhysicalInteractionLength(
                                                                 const G4Track& track,
                                                                 G4double   previousStepSize,
                                                                 G4ForceCondition* condition
                                                                 ) override;

    virtual inline G4VParticleChange* PostStepDoIt(
                                                   const G4Track& ,
                                                   const G4Step&
                                                   ) override;

    virtual inline G4double AtRestGetPhysicalInteractionLength(
                                                               const G4Track& ,
                                                               G4ForceCondition*
                                                               ) override;

    virtual inline G4VParticleChange* AtRestDoIt(
                                                 const G4Track& ,
                                                 const G4Step&
                                                 ) override;

    virtual inline G4double AlongStepGetPhysicalInteractionLength(
                                                                  const G4Track&,
                                                                  G4double  ,
                                                                  G4double  ,
                                                                  G4double& ,
                                                                  G4GPILSelection*
                                                                  ) override;

    virtual inline G4VParticleChange* AlongStepDoIt(
                                                    const G4Track& ,
                                                    const G4Step&
                                                    ) override;

  private:

    Mu2eSpecialCutsProcess(const Mu2eSpecialCutsProcess&) = delete;
    Mu2eSpecialCutsProcess(Mu2eSpecialCutsProcess&&) = delete;
    Mu2eSpecialCutsProcess& operator=(const Mu2eSpecialCutsProcess&) = delete;
    Mu2eSpecialCutsProcess& operator=(Mu2eSpecialCutsProcess&&) = delete;

    const Mu2eG4ResourceLimits& mu2elimits_;

  };
}
#endif
