#ifndef G4TNuclearCapture_h
#define G4TNuclearCapture_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Kevin Lynch, February 12 2010
// ----------------------------------------------------------------

#include "globals.hh"
#include "G4VRestDiscreteProcess.hh"
#include "G4ParticleChange.hh"

#include <vector>
// FIXME!  This can probably be derived from G4Decay, as it leverages
// almost ALL of the implementation.  Doesn't it?   .... actually, I
// think I should TEMPLATIZE this, and generate a
// G4MuAtomNuclearCapture and G4MuMoleculeNuclearCapture by typedef
template<class T> class G4TNuclearCapture : public G4VRestDiscreteProcess {
  
public:
  G4TNuclearCapture(G4String const& name = "TNuclearCapture");
  virtual ~G4TNuclearCapture();

  virtual G4double 
  PostStepGetPhysicalInteractionLength(const G4Track&, G4double,
				       G4ForceCondition*);
  virtual G4VParticleChange*
  PostStepDoIt(const G4Track& aTrack, const G4Step& aStep) {
    return CaptureIt(aTrack, aStep);
  }

  virtual G4double AtRestGetPhysicalInteractionLength(const G4Track& ,
						      G4ForceCondition*);
  virtual G4VParticleChange* 
  AtRestDoIt(const G4Track& aTrack, const G4Step&  aStep) {
    return CaptureIt(aTrack, aStep);
  }

  virtual G4bool IsApplicable(const G4ParticleDefinition&);

  virtual void StartTracking(G4Track*);
  virtual void EndTracking();

  void  SetVerboseLevel(G4int value){ verboseLevel = value; }
  G4int GetVerboseLevel() const { return verboseLevel; }

  G4double GetRemainderLifetime() const { return fRemainderLifeTime; }

protected:
  virtual G4double GetMeanFreePath(const G4Track& aTrack, G4double previousStepSize,
				   G4ForceCondition* condition);
  virtual G4double GetMeanLifeTime(const G4Track& aTrack,G4ForceCondition* condition);
  
  virtual G4VParticleChange* CaptureIt(const G4Track& aTrack, const G4Step& aStep);

  G4int verboseLevel;
  G4double fRemainderLifeTime;
  // FIXME: There's a pParticleChange in G4VProcess
  G4ParticleChange fParticleChange;

private:
  // disable these ... FIXME why?
  G4TNuclearCapture(G4TNuclearCapture<T> const&);
  G4TNuclearCapture& operator=(G4TNuclearCapture<T> const&);

};

#define G4TNUCLEARCAPTURE_ICC
#include "G4TNuclearCapture.icc"
#undef G4TNUCLEARCAPTURE_ICC

#endif
