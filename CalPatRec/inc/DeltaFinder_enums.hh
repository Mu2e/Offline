#ifndef DeltaFinder_enums_hh
#define DeltaFinder_enums_hh

#include "Offline/DataProducts/inc/StrawId.hh"

namespace mu2e {

  enum {
    kNStations      = StrawId::_nplanes/2,   // number of tracking stations
    kNFaces         = StrawId::_nfaces*2 ,   // N(faces) per station (4)
    kNPanelsPerFace = StrawId::_npanels/2    // = 3
  };

}
#endif
