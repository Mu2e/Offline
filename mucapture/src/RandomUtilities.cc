// ------------------------------------------------------------
//      GEANT 4 implementation file
//
//      History: first implementation, Kevin Lynch, March 18 2010
// ----------------------------------------------------------------

#include "RandomUtilities.hh"
#include "Randomize.hh"

G4ThreeVector GetRandomVec() {
  //
  // generate uniform vector
  //
  G4double const cost = 2.0 * G4UniformRand() - 1.0;
  G4double const sint = std::sqrt((1.0 - cost)*(1.0 + cost));
  G4double const Phi  = twopi * G4UniformRand();
  return G4ThreeVector(sint * std::cos(Phi), sint * std::sin(Phi), cost);
}

