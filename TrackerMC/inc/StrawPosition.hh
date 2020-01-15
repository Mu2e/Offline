#ifndef TrackerMC_StrawPosition_hh
#define TrackerMC_StrawPosition_hh
#include "Math/Vector3D.h"
namespace mu2e {
  namespace TrackerMC {
    typedef ROOT::Math::Cylindrical3D<Float_t> StrawPosition; // position WRT the straw.
  }
}
#endif
