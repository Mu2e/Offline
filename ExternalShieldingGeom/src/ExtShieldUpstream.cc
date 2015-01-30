#include "ExternalShieldingGeom/inc/ExtShieldUpstream.hh"

namespace mu2e {

  // genreflex persistency requires default ctr
  ExtShieldUpstream::ExtShieldUpstream() {}


  std::ostream& operator<<(std::ostream& os, const ExtShieldUpstream& ens) {
    os<<"ExtShieldUpstream("
      <<"material0="<<ens.getMaterialNames()[0]
//       <<", end vertices="<<ens.externalShieldOutline()
      <<", orientation0="<<ens.getOrientations()[0]
      <<" )";
    return os;
  }

} // namespace mu2e
