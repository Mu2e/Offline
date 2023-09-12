#ifndef STMGeom_InnerShielding_hh
#define STMGeom_InnerShielding_hh

// Germanium Detector Object
//
// Author: Haichuan Cao
// Sept 2023

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class InnerShielding {
  public:
    InnerShielding(bool build):
      _build(build)
    {}

    bool   build()                               const { return _build;  }


    InnerShielding() {}
  private:

    bool               _build;

  };

}

#endif/*STMGeom_InnerShielding_hh*/
