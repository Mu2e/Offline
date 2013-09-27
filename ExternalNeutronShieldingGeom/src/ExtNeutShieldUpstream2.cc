#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldUpstream2.hh"

namespace mu2e {

  // genreflex perstency requires default ctr
  ExtNeutShieldUpstream2::ExtNeutShieldUpstream2() {}


  std::ostream& operator<<(std::ostream& os, const ExtNeutShieldUpstream2& ens) {
    os<<"ExtNeutShieldUpstream2("
      <<"material="<<ens.materialName()
//       <<", end vertices="<<ens.externalShieldOutline()
      <<", length="<<ens.getLength()
      <<" )";
    return os;
  }

} // namespace mu2e
