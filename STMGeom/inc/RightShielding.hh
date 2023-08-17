#ifndef STMGeom_RightShielding_hh
#define STMGeom_RightShielding_hh

// Germanium Detector Object
//
// Author: Anthony Palladino
//

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class RightShielding {
  public:
    RightShielding(bool build):
      _build(build)
    {}

    bool   build()                               const { return _build;  }


    RightShielding() {}
  private:

    bool               _build;

  };

}

#endif/*STMGeom_RightShielding_hh*/
