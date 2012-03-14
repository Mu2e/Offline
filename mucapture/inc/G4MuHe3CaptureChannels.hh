#ifndef G4MuHe3CaptureChannels_hh
#define G4MuHe3CaptureChannels_hh 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Kevin Lynch, March 19 2010
// ----------------------------------------------------------------

#include "globals.hh"
#include "G4MuAtom.hh"
#include "G4VMuAtomCaptureKineticsChannel.hh"
#include "G4PhaseSpaceDecayChannel.hh"

class G4MuHe3ProtonChannel : public G4VMuAtomCaptureKineticsChannel {
public:
  G4MuHe3ProtonChannel(G4MuAtom const* p, G4double captureRate, G4int verboseLevel = 0);
  G4MuHe3ProtonChannel(G4MuAtom const* p, G4int verboseLevel = 0);
  G4DecayProducts* CaptureIt(G4DynamicParticle const*);
private:
  void CheckIsApplicable() const;
  G4PhaseSpaceDecayChannel psdc;
};

class G4MuHe3DeuteronChannel : public G4VMuAtomCaptureKineticsChannel {
public:
  G4MuHe3DeuteronChannel(G4MuAtom const* p, G4double captureRate, G4int verboseLevel = 0);
  G4MuHe3DeuteronChannel(G4MuAtom const* p, G4int verboseLevel = 0);
  G4DecayProducts* CaptureIt(G4DynamicParticle const*);
private:
  void CheckIsApplicable() const;
  G4PhaseSpaceDecayChannel psdc;
};

class G4MuHe3TritonChannel : public G4VMuAtomCaptureKineticsChannel {
public:
  G4MuHe3TritonChannel(G4MuAtom const* p, G4double captureRate, G4int verboseLevel = 0);
  G4MuHe3TritonChannel(G4MuAtom const* p, G4int verboseLevel = 0);
  G4DecayProducts* CaptureIt(G4DynamicParticle const*);
private:
  void CheckIsApplicable() const;
  G4PhaseSpaceDecayChannel psdc;
};


#endif
