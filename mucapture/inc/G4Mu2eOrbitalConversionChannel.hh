#ifndef G4Mu2eOrbitalConversionChannel_hh
#define G4Mu2eOrbitalConversionChannel_hh 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Kevin Lynch, March 23 2010
// ----------------------------------------------------------------

#include "globals.hh"
#include "G4MuAtom.hh"
#include "G4VMuAtomCaptureKineticsChannel.hh"
#include "G4PhaseSpaceDecayChannel.hh"

class G4Mu2eOrbitalConversionChannel : public G4VMuAtomCaptureKineticsChannel {
public:
  G4Mu2eOrbitalConversionChannel(G4MuAtom const* p, 
				 G4double captureRate, G4int verboseLevel = 0);
  G4DecayProducts* CaptureIt(G4DynamicParticle const*);
private:
  G4PhaseSpaceDecayChannel* GetPSDC(G4MuAtom const*);
  G4PhaseSpaceDecayChannel* psdc;
};

#endif
