#include "ProductionSolenoidGeom/inc/PSExternalShielding.hh"

namespace mu2e {

  // genreflex perstency requires default ctr
  PSExternalShielding::PSExternalShielding() {}


  std::ostream& operator<<(std::ostream& os, const PSExternalShielding& pse) {
    os<<"PSExternalShielding("
      <<"material="<<pse.materialName()
//       <<", end vertices="<<pse.externalShieldOutline()
      <<", length="<<pse.getLength()
      <<" )";
    return os;
  }

} // namespace mu2e
