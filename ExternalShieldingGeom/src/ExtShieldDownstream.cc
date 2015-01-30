#include "ExternalShieldingGeom/inc/ExtShieldDownstream.hh"

namespace mu2e {

  // genreflex persistency requires default ctr
  ExtShieldDownstream::ExtShieldDownstream() {}


  std::ostream& operator<<(std::ostream& os, const ExtShieldDownstream& ens) {
    os<<"ExtShieldDownstream("
      <<"Contains " << ens.getMaterialNames().size()
      <<" extrusions "
      <<" )";
    return os;
  }

} // namespace mu2e
