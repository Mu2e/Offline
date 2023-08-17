#ifndef STMGeom_BackShielding_hh
#define STMGeom_BackShielding_hh

// Germanium Detector Object
//
// Author: Anthony Palladino
//

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class BackShielding {
  public:
    BackShielding(bool build):
      _build(build)
    {}

    bool   build()                               const { return _build;  }


    BackShielding() {}
  private:

    bool               _build;

  };

}

#endif/*STMGeom_BackShielding_hh*/
