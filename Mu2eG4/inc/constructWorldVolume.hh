#ifndef Mu2eG4_constructWorldVolume_hh
#define Mu2eG4_constructWorldVolume_hh
//
// Free function to construct World Mother Volume
//
// $Id: constructWorldVolume.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
//
// Original author KLG
//

// Mu2e includes.
#include "G4Helper/inc/VolumeInfo.hh"

namespace mu2e {

  class SimpleConfig;

  VolumeInfo constructWorldVolume(SimpleConfig const * const _config
                            );

}

#endif /* Mu2eG4_constructWorldVolume_hh */
