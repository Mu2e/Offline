#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldCavexRight.hh"

namespace mu2e {

  // genreflex perstency requires default ctr
  ExtNeutShieldCavexRight::ExtNeutShieldCavexRight() {}


  std::ostream& operator<<(std::ostream& os, const ExtNeutShieldCavexRight& ens) {
    os<<"ExtNeutShieldCavexRight("
      <<"material="<<ens.materialName()
//       <<", end vertices="<<ens.externalShieldOutline()
      <<", length="<<ens.getLength()
      <<" )";
    return os;
  }

} // namespace mu2e
