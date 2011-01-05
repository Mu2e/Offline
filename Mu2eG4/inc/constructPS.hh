#ifndef constructPS_HH
#define constructPS_HH
//
// Free function to create  Production Solenoid and Production Target
//
// $Id: constructPS.hh,v 1.1 2011/01/05 21:04:31 genser Exp $
// $Author: genser $
// $Date: 2011/01/05 21:04:31 $
//
// Original author KLG
//

// G4 includes
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  void constructPS(VolumeInfo   const &       parent, 
                   SimpleConfig const * const _config,
                   G4ThreeVector&             _primaryProtonGunOrigin,
                   G4RotationMatrix&          _primaryProtonGunRotation
                   );

}

#endif
