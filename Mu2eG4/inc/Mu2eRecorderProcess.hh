#ifndef Mu2eG4_Mu2eRecorderProcess_h
#define Mu2eG4_Mu2eRecorderProcess_h 1

#include "G4VContinuousProcess.hh"

namespace mu2e{
  class Mu2eRecorderProcess : public G4VContinuousProcess
  {
  public:

    explicit Mu2eRecorderProcess(const G4String& processName ="Mu2eRecorderProcess" );

    virtual ~Mu2eRecorderProcess();

    virtual G4double GetContinuousStepLimit(const G4Track&, // track,
                                            G4double, // previousStepSize
                                            G4double, // currentMinimumStep
                                            G4double& // currentSafety
                                            );

    virtual G4VParticleChange* AlongStepDoIt(
                                             const G4Track& ,
                                             const G4Step&
                                             );

  private:

    // hide assignment operator as private
    Mu2eRecorderProcess& operator=(const Mu2eRecorderProcess&){return *this;};

  };
}
#endif

