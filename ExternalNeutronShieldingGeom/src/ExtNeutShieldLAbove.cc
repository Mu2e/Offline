#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldLAbove.hh"

namespace mu2e {

  // genreflex perstency requires default ctr
  ExtNeutShieldLAbove::ExtNeutShieldLAbove() {}


  std::ostream& operator<<(std::ostream& os, const ExtNeutShieldLAbove& ens) {
    os<<"ExtNeutShieldLAbove("
      <<"material="<<ens.materialName()
//       <<", end vertices="<<ens.externalShieldOutline()
      <<", length="<<ens.getLength()
      <<" )";
    return os;
  }

} // namespace mu2e
