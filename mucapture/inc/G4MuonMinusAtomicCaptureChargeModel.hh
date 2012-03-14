#ifndef G4MuonMinusAtomicCaptureChargeModel_h
#define G4MuonMinusAtomicCaptureChargeModel_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Kevin Lynch, February 22 2010
// ----------------------------------------------------------------

#include "globals.hh"
#include "G4Step.hh"

class G4MuonMinusAtomicCaptureChargeModel {

public:
  G4MuonMinusAtomicCaptureChargeModel(){}
  virtual ~G4MuonMinusAtomicCaptureChargeModel(){}

  // This should set the net MuAtom charge in units of eplus
  virtual G4int GetCharge(G4Step const&) const { return 0; }

protected:
};

#endif
