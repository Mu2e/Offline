#ifndef G4MuPStateModel_h
#define G4MuPStateModel_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Kevin Lynch, February 12 2010
// ----------------------------------------------------------------

#include "globals.hh"
#include "G4MuonMinusAtomicCaptureStateModel.hh"

class G4MuPStateModel : public G4MuonMinusAtomicCaptureStateModel {

public:
  G4MuPStateModel(){}
  virtual ~G4MuPStateModel(){}

  virtual G4int GetSpinState(G4Step const&) const;

protected:
};

#endif
