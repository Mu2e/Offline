#ifndef Mu2eG4_constructPS_hh
#define Mu2eG4_constructPS_hh
//
// Free function to create  Production Solenoid and Production Target
//
// $Id: constructPS.hh,v 1.5 2012/03/05 19:38:17 genser Exp $
// $Author: genser $
// $Date: 2012/03/05 19:38:17 $
//
// Original author KLG
//

// G4 includes
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  void constructPS(VolumeInfo const & parent, SimpleConfig const & _config);

}

#endif /* Mu2eG4_constructPS_hh */
