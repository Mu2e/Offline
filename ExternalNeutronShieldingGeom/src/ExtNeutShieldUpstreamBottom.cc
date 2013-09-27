#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldUpstreamBottom.hh"

namespace mu2e {

  // genreflex perstency requires default ctr
  ExtNeutShieldUpstreamBottom::ExtNeutShieldUpstreamBottom() {}


  std::ostream& operator<<(std::ostream& os, const ExtNeutShieldUpstreamBottom& ens) {
    os<<"ExtNeutShieldUpstreamBottom("
      <<"material="<<ens.materialName()
//       <<", end vertices="<<ens.externalShieldOutline()
      <<", length="<<ens.getLength()
      <<" )";
    return os;
  }

} // namespace mu2e
