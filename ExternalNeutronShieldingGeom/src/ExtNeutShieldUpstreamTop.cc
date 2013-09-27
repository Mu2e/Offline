#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldUpstreamTop.hh"

namespace mu2e {

  // genreflex perstency requires default ctr
  ExtNeutShieldUpstreamTop::ExtNeutShieldUpstreamTop() {}


  std::ostream& operator<<(std::ostream& os, const ExtNeutShieldUpstreamTop& ens) {
    os<<"ExtNeutShieldUpstreamTop("
      <<"material="<<ens.materialName()
//       <<", end vertices="<<ens.externalShieldOutline()
      <<", length="<<ens.getLength()
      <<" )";
    return os;
  }

} // namespace mu2e
