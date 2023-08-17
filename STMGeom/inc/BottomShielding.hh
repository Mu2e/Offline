#ifndef STMGeom_BottomShielding_hh
#define STMGeom_BottomShielding_hh

// Germanium Detector Object
//
// Author: Anthony Palladino
//

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class BottomShielding {
  public:
    BottomShielding(bool build):
      _build(build)
    {}

    bool   build()                               const { return _build;  }


    BottomShielding() {}
  private:

    bool               _build;

  };

}

#endif/*STMGeom_BottomShielding_hh*/
