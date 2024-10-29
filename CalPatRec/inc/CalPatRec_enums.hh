#ifndef CalPatRec_inc__enums_hh
#define CalPatRec_inc__enums_hh

#include "Offline/DataProducts/inc/StrawId.hh"

namespace mu2e {

  enum {
    kNStations      = StrawId::_nplanes/2,   // number of tracking stations
    kNFaces         = StrawId::_nfaces*2 ,   // N(faces) per station (4)
    kNPanelsPerFace = StrawId::_npanels/2    // = 3
  };

}
#endif
