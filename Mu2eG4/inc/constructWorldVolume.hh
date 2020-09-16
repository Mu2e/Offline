#ifndef Mu2eG4_constructWorldVolume_hh
#define Mu2eG4_constructWorldVolume_hh
//
// Free function to construct World Mother Volume
//
//
// Original author KLG
//

// Mu2e includes.
#include "G4Helper/inc/VolumeInfo.hh"

namespace mu2e {

  class SimpleConfig;
  VolumeInfo constructWorldVolume(const SimpleConfig& config);

}

#endif /* Mu2eG4_constructWorldVolume_hh */
