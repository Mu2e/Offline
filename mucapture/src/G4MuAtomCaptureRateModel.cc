// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Kevin Lynch, February 12 2010
// ----------------------------------------------------------------

#include "G4MuAtomCaptureRateModel.hh"
#include "G4StopElementSelector.hh"

G4double G4MuAtomCaptureRateModel::GetCaptureRate(G4int Z, G4int A, G4int /*iSpin*/) const {

  static G4StopElementSelector ses;

  return ses.GetMuonCaptureRate(Z,A);
}
