#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldCryoBoxes.hh"

namespace mu2e {

  // genreflex perstency requires default ctr
  ExtNeutShieldCryoBoxes::ExtNeutShieldCryoBoxes() {}


  std::ostream& operator<<(std::ostream& os, const ExtNeutShieldCryoBoxes& ens) {
    os<<"ExtNeutShieldCryoBoxes("
      <<"Number of boxes defined is:  "<<ens.materialNames().size()
      <<" )";
    return os;
  }

} // namespace mu2e
