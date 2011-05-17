#ifndef Mu2eG4_constructPS_hh
#define Mu2eG4_constructPS_hh
//
// Free function to create  Production Solenoid and Production Target
//
// $Id: constructPS.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
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
