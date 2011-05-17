#ifndef Mu2eG4_constructDirt_hh
#define Mu2eG4_constructDirt_hh
//
// Free function to create the earthen overburden.
//
// $Id: constructDirt.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
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
