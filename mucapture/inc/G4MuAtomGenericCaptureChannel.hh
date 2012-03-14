#ifndef G4MuAtomGenericCaptureChannel_h
#define G4MuAtomGenericCaptureChannel_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Kevin Lynch, March 19 2010
// ----------------------------------------------------------------

#include "globals.hh"
#include "G4MuAtom.hh"
#include "G4VMuAtomCaptureKineticsChannel.hh"

// FIXME ... does this use the capture rate model properly?
class G4MuAtomGenericCaptureChannel : public G4VMuAtomCaptureKineticsChannel {
public:
  G4MuAtomGenericCaptureChannel(G4MuAtom const* p, G4int verboseLevel = 0);
  G4MuAtomGenericCaptureChannel(G4MuAtom const* p, G4double rate, G4int verboseLevel = 0); 
  G4bool IsCloneable() const { return true; }
  G4MuAtomGenericCaptureChannel* Clone(G4MuAtom const*);
  G4DecayProducts* CaptureIt(G4DynamicParticle const*);
};

#endif
