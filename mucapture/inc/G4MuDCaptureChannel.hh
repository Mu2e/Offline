#ifndef G4MuDCaptureChannel_hh
#define G4MuDCaptureChannel_hh 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Kevin Lynch, March 18 2010
// ----------------------------------------------------------------

#include "globals.hh"
#include "G4VMuAtomCaptureKineticsChannel.hh"
#include "G4PhaseSpaceDecayChannel.hh"

class G4MuDCaptureChannel : public G4VMuAtomCaptureKineticsChannel {
public:
  G4MuDCaptureChannel(G4MuAtom const* p, G4double captureRate, G4int verboseLevel = 0);
  G4MuDCaptureChannel(G4MuAtom const* p, G4int verboseLevel = 0);
  G4DecayProducts* CaptureIt(G4DynamicParticle const*);
private:
  void CheckIsApplicable() const;
  G4PhaseSpaceDecayChannel psdc;
};

#endif
