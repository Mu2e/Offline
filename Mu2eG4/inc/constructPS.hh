#ifndef Mu2eG4_constructPS_hh
#define Mu2eG4_constructPS_hh
//
// Free function to create  Production Solenoid and Production Target
//
//
// Original author KLG
//

//Mu2e includes
#include "Mu2eG4Helper/inc/VolumeInfo.hh"

class G4LogicalVolume;

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  G4LogicalVolume* constructPS(VolumeInfo const & parent, SimpleConfig const & _config);

}

#endif /* Mu2eG4_constructPS_hh */
