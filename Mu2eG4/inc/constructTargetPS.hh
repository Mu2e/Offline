#ifndef Mu2eG4_constructTargetPS_hh
#define Mu2eG4_constructTargetPS_hh
//
// Free function to create Production Target
//
//
// Original author Giovanni F. Tassielli
//

// G4 includes
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4RotationMatrix.hh"

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  void constructTargetPS(VolumeInfo const & parent, SimpleConfig const & _config);

}

#endif /* Mu2eG4_constructTargetPS_hh */
