#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldLCeiling.hh"

namespace mu2e {

  // genreflex perstency requires default ctr
  ExtNeutShieldLCeiling::ExtNeutShieldLCeiling() {}


  std::ostream& operator<<(std::ostream& os, const ExtNeutShieldLCeiling& ens) {
    os<<"ExtNeutShieldLCeiling("
      <<"material="<<ens.materialName()
//       <<", end vertices="<<ens.externalShieldOutline()
      <<", length="<<ens.getLength()
      <<" )";
    return os;
  }

} // namespace mu2e
