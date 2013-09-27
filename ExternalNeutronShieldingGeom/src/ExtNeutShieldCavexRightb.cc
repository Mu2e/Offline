#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldCavexRightb.hh"

namespace mu2e {

  // genreflex perstency requires default ctr
  ExtNeutShieldCavexRightb::ExtNeutShieldCavexRightb() {}


  std::ostream& operator<<(std::ostream& os, const ExtNeutShieldCavexRightb& ens) {
    os<<"ExtNeutShieldCavexRightb("
      <<"material="<<ens.materialName()
//       <<", end vertices="<<ens.externalShieldOutline()
      <<", length="<<ens.getLength()
      <<" )";
    return os;
  }

} // namespace mu2e
