#ifndef Mu2eG4_constructPS_hh
#define Mu2eG4_constructPS_hh
//
// Free function to create  Production Solenoid and Production Target
//
// $Id: constructPS.hh,v 1.3 2011/05/18 02:27:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:17 $
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

#endif /* Mu2eG4_constructPS_hh */
