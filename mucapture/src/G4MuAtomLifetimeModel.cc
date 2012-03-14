// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Kevin Lynch, February 22 2010
// ----------------------------------------------------------------

#include "G4MuAtomLifetimeModel.hh"
#include "G4MuonPlus.hh"

// Updated from G4StopElementSelector::GetMuonDecayRate
// Decay time on K-shell 
// N.C.Mukhopadhyay Phys. Rep. 30 (1977) 1.
G4double G4MuAtomLifetimeModel::GetLifetime(G4int Z, G4int /*A*/) const {
  G4double const alphaZ = Z/137.;
  G4double lambda = 1.0 - 2.5*alphaZ*alphaZ ;
  if( 0.5 > lambda ) lambda = 0.5;
  return G4MuonPlus::Definition()->GetPDGLifeTime()/lambda;
}

