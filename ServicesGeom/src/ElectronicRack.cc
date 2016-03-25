#include "ServicesGeom/inc/ElectronicRack.hh"

namespace mu2e {

  // genreflex persistency requires default ctr
  ElectronicRack::ElectronicRack() {}


  std::ostream& operator<<(std::ostream& os, const ElectronicRack& ens) {
    os<<"ElectronicRack("
      << "Contains " << ens.getMaterialNames().size()
      <<" Rack enclosures "
      <<" )";
    return os;
  }

} // namespace mu2e
