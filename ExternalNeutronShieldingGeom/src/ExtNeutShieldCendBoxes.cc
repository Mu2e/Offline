#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldCendBoxes.hh"

namespace mu2e {

  // genreflex perstency requires default ctr
  ExtNeutShieldCendBoxes::ExtNeutShieldCendBoxes() {}


  std::ostream& operator<<(std::ostream& os, const ExtNeutShieldCendBoxes& ens) {
    os<<"ExtNeutShieldCendBoxes("
      <<"Number of boxes defined is:  "<<ens.materialNames().size()
      <<" )";
    return os;
  }

} // namespace mu2e
