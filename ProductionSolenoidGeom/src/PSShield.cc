#include "ProductionSolenoidGeom/inc/PSShield.hh"

namespace mu2e {

  std::ostream& operator<<(std::ostream& os, const PSShield& shield) {
    return os<<"PSShield("
             <<shield.bulk()
      //<<", originInMu2e="<<shield.originInMu2e()
             <<", cutoutRefPoint="<<shield.cutoutRefPoint()
             <<", cutoutR="<<shield.cutoutR()
             <<", cutoutHalfLength="<<shield.cutoutHalfLength()
             <<", cutoutRotY="<<shield.cutoutRotY()
             <<")";
  }

} // namespace mu2e
