#include "RecoDataProducts/inc/ExtMonFNALRawCluster.hh"

namespace mu2e {

  std::ostream& operator<<(std::ostream& os, const ExtMonFNALRawCluster& c) {
    os<<"ExtMonFNALRawCluster(nhits = "<<c.hits().size()
      <<", hits = { ";
    for(ExtMonFNALRawCluster::Hits::const_iterator i = c.hits().begin(); i != c.hits().end(); ++i) {
      os<<**i<<", ";
    }
    os<<" } )";
    return os;
  }

} // namespace mu2e
