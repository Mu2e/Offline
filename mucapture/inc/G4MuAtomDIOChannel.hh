#ifndef G4MuAtomDIOChannel_h
#define G4MuAtomDIOChannel_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Kevin Lynch, February 12 2010
// ----------------------------------------------------------------

#include "G4VMuAtomDecayChannel.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class G4MuAtomDIOChannel : public G4VMuAtomDecayChannel {
public:

  G4MuAtomDIOChannel(G4String const& parentName,
		     G4double BR);

  virtual G4MuAtomDIOChannel* Clone();
  virtual ~G4MuAtomDIOChannel();
  virtual G4DecayProducts* DecayIt(G4double);

};

#endif
