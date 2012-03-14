#ifndef G4MuonMinusAtomicCaptureStateModel_h
#define G4MuonMinusAtomicCaptureStateModel_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Kevin Lynch, February 12 2010
// ----------------------------------------------------------------

#include "globals.hh"
#include "G4Step.hh"

class G4MuonMinusAtomicCaptureStateModel {

public:
  G4MuonMinusAtomicCaptureStateModel(){}
  virtual ~G4MuonMinusAtomicCaptureStateModel(){}

  virtual G4int GetSpinState(G4Step const&) const { return 0; }

protected:
};

#endif
