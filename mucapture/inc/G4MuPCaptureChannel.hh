#ifndef G4MuPCaptureChannel_hh
#define G4MuPCaptureChannel_hh 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Kevin Lynch, March 12 2010
// ----------------------------------------------------------------

#include "globals.hh"
#include "G4MuAtom.hh"
#include "G4VMuAtomCaptureKineticsChannel.hh"

class G4MuPCaptureChannel : public G4VMuAtomCaptureKineticsChannel {
public:
  G4MuPCaptureChannel(G4MuAtom const* p, G4double captureRate, G4int verboseLevel = 0);
  G4MuPCaptureChannel(G4MuAtom const* p, G4int verboseLevel = 0);
  G4DecayProducts* CaptureIt(G4DynamicParticle const*);
private:
  void CheckIsApplicable() const;
};

#endif
