#ifndef GeometryService_PTMonMaker_hh
#define GeometryService_PTMonMaker_hh

#include <memory>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "ConfigTools/inc/SimpleConfig.hh"
#include "PTMonGeom/inc/PTMonPWC.hh"
#include "PTMonGeom/inc/PTMon.hh"
//
// construct and return a PTMon
//
// original author Helenka Casler
//

namespace mu2e {

  class PTMonMaker {

  public:
    static std::unique_ptr<PTMon> make(SimpleConfig const& config);
  };

} // namespace mu2e

#endif