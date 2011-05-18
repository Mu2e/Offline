#ifndef Mu2eG4_constructDirt_hh
#define Mu2eG4_constructDirt_hh
//
// Free function to create the earthen overburden.
//
// $Id: constructDirt.hh,v 1.3 2011/05/18 02:27:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:17 $
//
// Original author KLG
//

// Mu2e includes.
#include "G4Helper/inc/VolumeInfo.hh"

namespace mu2e {

  class SimpleConfig;

  VolumeInfo constructDirt(VolumeInfo   const & parent,
                           SimpleConfig const * const _config
                           );

}

#endif /* Mu2eG4_constructDirt_hh */
