#ifndef STMGeom_LeftShielding_hh
#define STMGeom_LeftShielding_hh

// Germanium Detector Object
//
// Author: Anthony Palladino
//

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class LeftShielding {
  public:
    LeftShielding(bool build):
      _build(build)
    {}

    bool   build()                               const { return _build;  }


    LeftShielding() {}
  private:

    bool               _build;

  };

}

#endif/*STMGeom_LeftShielding_hh*/
