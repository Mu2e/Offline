#ifndef Mu2eG4_constructTargetPS_hh
#define Mu2eG4_constructTargetPS_hh
//
// Free function to create Production Target
//
// $Id: constructTargetPS.hh,v 1.1 2013/09/06 19:39:13 tassiell Exp $
// $Author: tassiell $
// $Date: 2013/09/06 19:39:13 $
//
// Original author Giovanni F. Tassielli
//

// G4 includes
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  void constructTargetPS(VolumeInfo const & parent, SimpleConfig const & _config);

}

#endif /* Mu2eG4_constructTargetPS_hh */
