// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Kevin Lynch, February 12 2010
// ----------------------------------------------------------------

#include "G4MuAtomCaptureRateModel.hh"
#include "G4MuonMinusBoundDecay.hh"

G4double G4MuAtomCaptureRateModel::GetCaptureRate(G4int Z, G4int A, G4int /*iSpin*/) const {
  return G4MuonMinusBoundDecay::GetMuonCaptureRate(Z,A);
}
