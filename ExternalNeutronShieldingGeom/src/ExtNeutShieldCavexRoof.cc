#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldCavexRoof.hh"

namespace mu2e {

  // genreflex perstency requires default ctr
  ExtNeutShieldCavexRoof::ExtNeutShieldCavexRoof() {}


  std::ostream& operator<<(std::ostream& os, const ExtNeutShieldCavexRoof& ens) {
    os<<"ExtNeutShieldCavexRoof("
      <<"material="<<ens.materialName()
//       <<", end vertices="<<ens.externalShieldOutline()
      <<", length="<<ens.getLength()
      <<" )";
    return os;
  }

} // namespace mu2e
