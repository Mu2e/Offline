#include "RecoDataProducts/inc/ExtMonFNALRecoCluster.hh"

namespace mu2e {

  ExtMonFNALRecoCluster::ExtMonFNALRecoCluster(const art::Ptr<ExtMonFNALRawCluster>& raw,
                                               unsigned int plane,
                                               const CLHEP::Hep3Vector& position,
                                               int xWidth,
                                               int yWidth,
                                               int clock)
    : raw_(raw)
    , position_(position)
    , plane_(plane)
    , xWidth_(xWidth)
    , yWidth_(yWidth)
    , clock_(clock)
  {}

  std::ostream& operator<<(std::ostream& os, const ExtMonFNALRecoCluster& c) {
    os<<"ExtMonFNALRecoCluster(plane="<<c.plane()
      <<", pos="<<c.position()
      <<", xw="<<c.xWidth()
      <<", yw="<<c.yWidth()
      <<", clock="<<c.clock()
      <<" )";
      ;
    return os;
  }

} // namespace mu2e
