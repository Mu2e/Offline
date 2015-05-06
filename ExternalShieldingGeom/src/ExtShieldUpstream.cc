#include "ExternalShieldingGeom/inc/ExtShieldUpstream.hh"

namespace mu2e {

  // genreflex persistency requires default ctr
  ExtShieldUpstream::ExtShieldUpstream() {}


  std::ostream& operator<<(std::ostream& os, const ExtShieldUpstream& ens) {
    os<<"ExtShieldUpstream("
      << "Contains " << ens.getMaterialNames().size()
      <<" boxes "
      <<" )";
    return os;
  }

} // namespace mu2e
