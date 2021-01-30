#ifndef Mu2eG4_Mu2eGammaDaughterProcess_h
#define Mu2eG4_Mu2eGammaDaughterProcess_h

#include "G4VDiscreteProcess.hh"

namespace mu2e{
  class Mu2eGammaDaughterProcess : public G4VDiscreteProcess
  {
  public:

    Mu2eGammaDaughterProcess(const G4String& processName ="Mu2eGammaDaughterProcess" );

    virtual ~Mu2eGammaDaughterProcess();

    virtual G4bool IsApplicable(const G4ParticleDefinition&);


    virtual G4double PostStepGetPhysicalInteractionLength( const G4Track& track,
                                                           G4double previousStepSize,
                                                           G4ForceCondition* condition);

    virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

    virtual G4double GetMeanFreePath(const G4Track&, G4double,G4ForceCondition*);

    void SetMinDaughterEnergy(G4double energy) { minDaughterEnergy_ = energy; }
    void SetKillAfterConvert(G4bool killAfterConvert) { killAfterConvert_ = killAfterConvert; }
    void SetVerboseLevel(G4int verbose) { verbose_ = verbose; }

  private:

    // hide assignment operator as private
    Mu2eGammaDaughterProcess& operator=(const Mu2eGammaDaughterProcess&){return *this;};

    G4double photonEnergy_;
    G4double minDaughterEnergy_;
    G4bool killAfterConvert_;
    G4int accepted_;
    G4int verbose_;
  };
}
#endif
