#ifndef Mu2eG4_constructHall_hh
#define Mu2eG4_constructHall_hh
//
// Free function to create the hall walls and hall interior inside the earthen overburden.
//
// $Id: constructHall.hh,v 1.4 2012/04/17 19:56:56 gandr Exp $
// $Author: gandr $
// $Date: 2012/04/17 19:56:56 $
//
// Original author KLG
//

// Mu2e includes.
#include "G4Helper/inc/VolumeInfo.hh"

namespace mu2e {

  class SimpleConfig;

  VolumeInfo constructHall(const VolumeInfo& worldInfo, const SimpleConfig& config);

}

#endif /* Mu2eG4_constructHall_hh */
