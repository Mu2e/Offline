#include "RecoDataProducts/inc/ExtMonFNALTrkParam.hh"

namespace mu2e {

  std::ostream& operator<<(std::ostream& os, const ExtMonFNALTrkParam& p) {
    os<<"ExtMonFNALTrkParam(z0="<<p.z0()
      <<", x="<<p.posx()
      <<", sx="<<p.slopex()
      <<", y="<<p.posy()
      <<", sy="<<p.slopey()
      <<", rinv="<<p.rinv()
      <<" )";
      ;
    return os;
  }

} // namespace mu2e
