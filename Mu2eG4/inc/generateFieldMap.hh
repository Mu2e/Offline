#ifndef Mu2eG4_generateFieldMap_hh
#define Mu2eG4_generateFieldMap_hh
//
// Free function to get field, used by G4, at particular points in space.
// Useful for different tests.
//
// Original author Ivan Logashenko
//

// Geant4 includes
#include "Geant4/G4ThreeVector.hh"

namespace mu2e {

  void generateFieldMap(const G4ThreeVector &, int);

}

#endif /* Mu2eG4_generateFieldMap_hh */
