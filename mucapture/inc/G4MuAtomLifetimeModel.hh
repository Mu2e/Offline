#ifndef G4MuAtomLifetimeModel_h
#define G4MuAtomLifetimeModel_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Kevin Lynch, February 12 2010
// ----------------------------------------------------------------

#include "globals.hh"

class G4MuAtomLifetimeModel {
public:
  G4MuAtomLifetimeModel(){}
  virtual ~G4MuAtomLifetimeModel(){}

  virtual G4double GetLifetime(G4int Z, G4int A) const;

};

#endif
