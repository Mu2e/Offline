// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Kevin Lynch, February 12 2010
// ----------------------------------------------------------------

#include "G4MuPStateModel.hh"
#include "Randomize.hh"

G4int G4MuPStateModel::GetSpinState(G4Step const&) const {
  G4double const val = G4UniformRand();
  G4int ret;
  if( val < 3./4. )
    ret = 2;
  else
    ret = 0;
  return ret;
}
