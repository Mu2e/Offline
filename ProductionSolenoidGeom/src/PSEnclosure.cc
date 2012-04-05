#include "ProductionSolenoidGeom/inc/PSEnclosure.hh"

namespace mu2e {

  // genreflex perstency requires default ctr
  PSEnclosure::PSEnclosure()
    : shell_(Tube::UNINITIALIZED),
      vacuum_(Tube::UNINITIALIZED),
      endPlate_(Tube::UNINITIALIZED)
  {}

} // namespace mu2e
