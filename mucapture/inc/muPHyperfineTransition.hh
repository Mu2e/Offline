#ifndef muPHyperfineTransition_hh
#define muPHyperfineTransition_hh 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Kevin Lynch, March 12 2010
// ----------------------------------------------------------------

#include "globals.hh"
#include "G4VMuAtomCaptureKineticsChannel.hh"

class muPHyperfineStoT : public G4VMuAtomCaptureKineticsChannel {
public:
  muPHyperfineStoT(G4MuAtom const* p, G4double transitionRate, G4int verboseLevel = 0);
  G4DecayProducts* CaptureIt(G4DynamicParticle const*);
private:
  void CheckIsApplicable() const;
};

class muPHyperfineTtoS : public G4VMuAtomCaptureKineticsChannel {
public:
  muPHyperfineTtoS(G4MuAtom const* p, G4double transitionRate, G4int verboseLevel = 0);
  G4DecayProducts* CaptureIt(G4DynamicParticle const*);
private:
  void CheckIsApplicable() const;
};



#endif
