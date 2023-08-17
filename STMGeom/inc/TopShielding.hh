#ifndef STMGeom_TopShielding_hh
#define STMGeom_TopShielding_hh

// Germanium Detector Object
//
// Author: Anthony Palladino
//

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class TopShielding {
  public:
    TopShielding(bool build):
      _build(build)
    {}

    bool   build()                               const { return _build;  }


    TopShielding() {}
  private:

    bool               _build;

  };

}

#endif/*STMGeom_TopShielding_hh*/
