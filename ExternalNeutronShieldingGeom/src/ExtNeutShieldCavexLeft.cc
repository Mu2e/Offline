#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldCavexLeft.hh"
// This code will be superceded in a package called ExternalShieldingGeom
namespace mu2e {

  // genreflex perstency requires default ctr
  ExtNeutShieldCavexLeft::ExtNeutShieldCavexLeft() {}


  std::ostream& operator<<(std::ostream& os, const ExtNeutShieldCavexLeft& ens) {
    os<<"ExtNeutShieldCavexLeft("
      <<"material="<<ens.materialName()
//       <<", end vertices="<<ens.externalShieldOutline()
      <<", length="<<ens.getLength()
      <<" )";
    return os;
  }

} // namespace mu2e
