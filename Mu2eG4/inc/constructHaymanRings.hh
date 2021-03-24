#ifndef Mu2eG4_constructHaymanRings_hh
#define Mu2eG4_constructHaymanRings_hh
//
// Free function to create strawman Production Target

// G4 includes
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4RotationMatrix.hh"
#include "TMath.h"
namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  void constructHaymanRings(VolumeInfo const & parent, SimpleConfig const & _config);

}

#endif /* Mu2eG4_constructHaymanRings_hh */
