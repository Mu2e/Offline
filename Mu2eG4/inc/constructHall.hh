#ifndef constructHall_HH
#define constructHall_HH
//
// Free function to create the hall walls and hall interior inside the earthen overburden.
//
// $Id: constructHall.hh,v 1.1 2011/01/05 21:04:31 genser Exp $
// $Author: genser $
// $Date: 2011/01/05 21:04:31 $
//
// Original author KLG
//

// Mu2e includes.
#include "G4Helper/inc/VolumeInfo.hh"

namespace mu2e {

  class SimpleConfig;

  VolumeInfo constructHall(VolumeInfo   const & parent, 
                           SimpleConfig const * const _config
                           );

}

#endif
