#include "ProductionSolenoidGeom/inc/PSShield.hh"

namespace mu2e {

  // genreflex perstency requires default ctr
  PSShield::PSShield()
    : bulk_(std::vector<double>(), std::vector<double>(), std::vector<double>(), CLHEP::Hep3Vector(), "")
    , cutoutRefPoint_(), cutoutR_(), cutoutHalfLength_(), cutoutRotY_()
  {}


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
