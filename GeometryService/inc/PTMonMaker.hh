#ifndef GeometryService_PTMonMaker_hh
#define GeometryService_PTMonMaker_hh

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"
//
// construct and return a PTMon
//
// original author Helenka Casler
//

namespace mu2e {

  class SimpleConfig;
  class PTMon;
  // TODO: class PTMon (or something)

  class PTMonMaker {

  public:
    static std::unique_ptr<PTMon> make(SimpleConfig const& config);
  };

} // namespace mu2e

#endif