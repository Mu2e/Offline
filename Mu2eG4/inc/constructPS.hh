#ifndef Mu2eG4_constructPS_hh
#define Mu2eG4_constructPS_hh
//
// Free function to create  Production Solenoid and Production Target
//
// $Id: constructPS.hh,v 1.4 2011/11/02 21:29:27 gandr Exp $
// $Author: gandr $
// $Date: 2011/11/02 21:29:27 $
//
// Original author KLG
//

// G4 includes
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  void constructPS(const VolumeInfo& parent, const SimpleConfig * const _config);

}

#endif /* Mu2eG4_constructPS_hh */
