#include "ExternalShieldingGeom/inc/Saddle.hh"

namespace mu2e {

  // genreflex persistency requires default ctr
  Saddle::Saddle() {}


  std::ostream& operator<<(std::ostream& os, const Saddle& ens) {
    os<<"Saddle("
      <<"Contains " << ens.getMaterialNames().size()
      <<" extrusions "
      <<" )";
    return os;
  }

} // namespace mu2e
