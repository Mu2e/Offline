#ifndef GeometryService_PTMMaker_hh
#define GeometryService_PTMMaker_hh

#include <memory>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/PTMGeom/inc/PTMPWC.hh"
#include "Offline/PTMGeom/inc/PTM.hh"
//
// construct and return a PTMon
//
// original author Helenka Casler
//

namespace mu2e {

  class PTMMaker {

  public:
    static std::unique_ptr<PTM> make(SimpleConfig const& config);
    static std::unique_ptr<PTM> makeWithBasicStand(SimpleConfig const& config);
    static std::unique_ptr<PTM> makeFloatingPWCs(SimpleConfig const& config);
  };

} // namespace mu2e

#endif
