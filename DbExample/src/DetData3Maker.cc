#include "DbExample/inc/DetData3Maker.hh"

namespace mu2e {

  DetData3::ptr_t DetData3Maker::fromConfig(SimpleConfig const& sc) {
    std::vector<double> vv;
    std::vector<double> vvd(3,-2.0);
    sc.getVectorDouble ("response",vv,vvd);
    return fromStatic(vv);
  }
  
  DetData3::ptr_t DetData3Maker::fromStatic(std::vector<double> const& vv) {
    std::array<double,3> aa;
    for(int i=0; i<3; i++) aa[i] = vv[i];
    return std::make_shared<DetData3>(aa);
  }
  
  DetData3::ptr_t DetData3Maker::fromDb(DetData1 const& d1, 
					DetData2 const& d2) {
    std::array<double,3> aa = d2.getData();
    for(int channel=0; channel<3; channel++) {
      aa[channel] *= d1.getData();
    }
    return std::make_shared<DetData3>(aa);
  }

}
