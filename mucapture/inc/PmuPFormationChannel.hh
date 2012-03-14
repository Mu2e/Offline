#ifndef PmuPFormationChannel_hh
#define PmuPFormationChannel_hh 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Kevin Lynch, March 12 2010
// ----------------------------------------------------------------

#include "globals.hh"
#include "G4MuAtom.hh"
#include "G4VMuAtomCaptureKineticsChannel.hh"

class PmuPFormationChannel : public G4VMuAtomCaptureKineticsChannel {
public:
  PmuPFormationChannel(G4MuAtom const* p, G4double captureRate, G4int verboseLevel = 0);
  //  PmuPFormationChannel(G4MuAtom const* p, G4int verboseLevel = 0);
  G4DecayProducts* CaptureIt(G4DynamicParticle const*);
private:
  void CheckIsApplicable() const;
};

#endif
