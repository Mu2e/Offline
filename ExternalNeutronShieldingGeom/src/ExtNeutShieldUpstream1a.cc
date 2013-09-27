#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldUpstream1a.hh"

namespace mu2e {

  // genreflex perstency requires default ctr
  ExtNeutShieldUpstream1a::ExtNeutShieldUpstream1a() {}


  std::ostream& operator<<(std::ostream& os, const ExtNeutShieldUpstream1a& ens) {
    os<<"ExtNeutShieldUpstream1a("
      <<"material="<<ens.materialName()
//       <<", end vertices="<<ens.externalShieldOutline()
      <<", length="<<ens.getLength()
      <<" )";
    return os;
  }

} // namespace mu2e
