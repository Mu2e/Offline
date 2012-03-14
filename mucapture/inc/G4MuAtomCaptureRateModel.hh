#ifndef G4MuAtomCaptureRateModel_h
#define G4MuAtomCaptureRateModel_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Kevin Lynch, February 12 2010
// ----------------------------------------------------------------

#include "globals.hh"

class G4MuAtomCaptureRateModel {
public:
  G4MuAtomCaptureRateModel(){}
  virtual ~G4MuAtomCaptureRateModel(){}

  virtual G4double GetCaptureRate(G4int Z, G4int A, G4int iSpin) const;

};

#endif
