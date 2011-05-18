#ifndef Mu2eG4_constructHall_hh
#define Mu2eG4_constructHall_hh
//
// Free function to create the hall walls and hall interior inside the earthen overburden.
//
// $Id: constructHall.hh,v 1.3 2011/05/18 02:27:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:17 $
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

#endif /* Mu2eG4_constructHall_hh */
