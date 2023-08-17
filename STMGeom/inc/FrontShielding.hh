#ifndef STMGeom_FrontShielding_hh
#define STMGeom_FrontShielding_hh

// Germanium Detector Object
//
// Author: Anthony Palladino
//

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class FrontShielding {
  public:
    FrontShielding(bool build):
      _build(build)
    {}

    bool   build()                               const { return _build;  }


    FrontShielding() {}
  private:

    bool               _build;

  };

}

#endif/*STMGeom_FrontShielding_hh*/
