#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldUpstream1b.hh"

namespace mu2e {

  // genreflex perstency requires default ctr
  ExtNeutShieldUpstream1b::ExtNeutShieldUpstream1b() {}


  std::ostream& operator<<(std::ostream& os, const ExtNeutShieldUpstream1b& ens) {
    os<<"ExtNeutShieldUpstream1b("
      <<"material="<<ens.materialName()
//       <<", end vertices="<<ens.externalShieldOutline()
      <<", length="<<ens.getLength()
      <<" )";
    return os;
  }

} // namespace mu2e
