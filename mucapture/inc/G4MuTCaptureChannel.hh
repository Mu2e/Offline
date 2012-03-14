#ifndef G4MuTCaptureChannel_hh
#define G4MuTCaptureChannel_hh 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Kevin Lynch, March 18 2010
// ----------------------------------------------------------------

#include "globals.hh"
#include "G4MuAtom.hh"
#include "G4VMuAtomCaptureKineticsChannel.hh"
#include "G4PhaseSpaceDecayChannel.hh"

class G4MuTCaptureChannel : public G4VMuAtomCaptureKineticsChannel {
public:
  G4MuTCaptureChannel(G4MuAtom const* p, G4double captureRate, G4int verboseLevel = 0);
  G4MuTCaptureChannel(G4MuAtom const* p, G4int verboseLevel = 0);
  G4DecayProducts* CaptureIt(G4DynamicParticle const*);
private:
  void CheckIsApplicable() const;
  G4PhaseSpaceDecayChannel psdc;
};

#endif
