#ifndef constructWorldVolume_HH
#define constructWorldVolume_HH
//
// Free function to construct World Mother Volume
//
// $Id: constructWorldVolume.hh,v 1.1 2011/01/05 21:04:31 genser Exp $
// $Author: genser $
// $Date: 2011/01/05 21:04:31 $
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

#endif
