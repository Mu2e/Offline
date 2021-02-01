#ifndef Mu2eG4_Mu2eGammaDaughterCut_h
#define Mu2eG4_Mu2eGammaDaughterCut_h

#include "G4VDiscreteProcess.hh"

namespace mu2e{
  class Mu2eGammaDaughterCut : public G4VDiscreteProcess
  {
  public:

    Mu2eGammaDaughterCut(const G4double minDaughterEnergy, const G4bool killAfterConvert = false,
                         const G4int verbose = 0,
                         const G4String& processName ="Mu2eGammaDaughterCut" );

    virtual ~Mu2eGammaDaughterCut();

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
    Mu2eGammaDaughterCut& operator=(const Mu2eGammaDaughterCut&){return *this;};

    G4double minDaughterEnergy_;
    G4bool killAfterConvert_;
    G4int verbose_;
    G4int accepted_;
    G4double photonEnergy_;
  };
}
#endif
