#ifndef G4MuHe4CaptureChannels_hh
#define G4MuHe4CaptureChannels_hh 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Kevin Lynch, March 19 2010
// ----------------------------------------------------------------

#include "globals.hh"
#include "G4MuAtom.hh"
#include "G4VMuAtomCaptureKineticsChannel.hh"
#include "G4PhaseSpaceDecayChannel.hh"

class G4MuHe4ProtonChannel : public G4VMuAtomCaptureKineticsChannel {
public:
  G4MuHe4ProtonChannel(G4MuAtom const* p, G4double captureRate, G4int verboseLevel = 0);
  G4MuHe4ProtonChannel(G4MuAtom const* p, G4int verboseLevel = 0);
  G4DecayProducts* CaptureIt(G4DynamicParticle const*);
private:
  void CheckIsApplicable() const;
  G4PhaseSpaceDecayChannel psdc;
};

class G4MuHe4DeuteronChannel : public G4VMuAtomCaptureKineticsChannel {
public:
  G4MuHe4DeuteronChannel(G4MuAtom const* p, G4double captureRate, G4int verboseLevel = 0);
  G4MuHe4DeuteronChannel(G4MuAtom const* p, G4int verboseLevel = 0);
  G4DecayProducts* CaptureIt(G4DynamicParticle const*);
private:
  void CheckIsApplicable() const;
  G4PhaseSpaceDecayChannel psdc;
};

class G4MuHe4TritonChannel : public G4VMuAtomCaptureKineticsChannel {
public:
  G4MuHe4TritonChannel(G4MuAtom const* p, G4double captureRate, G4int verboseLevel = 0);
  G4MuHe4TritonChannel(G4MuAtom const* p, G4int verboseLevel = 0);
  G4DecayProducts* CaptureIt(G4DynamicParticle const*);
private:
  void CheckIsApplicable() const;
  G4PhaseSpaceDecayChannel psdc;
};


#endif
