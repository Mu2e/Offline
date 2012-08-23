#include "RecoDataProducts/inc/ExtMonFNALRawHit.hh"

namespace mu2e {

  std::ostream& operator<<(std::ostream& os, const ExtMonFNALRawHit& hit) {
    return os<<"ExtMonFNALRawHit(pix="<<hit.pixelId()
             <<", clock = "<<hit.clock()
             <<", tot = "<<hit.tot()
             <<" )";
  }

} // namespace mu2e
