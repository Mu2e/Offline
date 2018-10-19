#include "DbExample/inc/DetData2Maker.hh"

namespace mu2e {

  DetData2::ptr_t DetData2Maker::fromConfig(SimpleConfig const& sc) {
    std::vector<double> vv;
    std::vector<double> vvd(3,-1.0);
    sc.getVectorDouble ("dtoe",vv,vvd);
    return fromStatic(vv);
  }
  
  DetData2::ptr_t DetData2Maker::fromStatic(std::vector<double> const& vv) {
    std::array<double,3> aa;
    for(int i=0; i<3; i++) aa[i] = vv[i];
    return std::make_shared<DetData2>(aa);
  }
  
  DetData2::ptr_t DetData2Maker::fromDb(TstCalib1 const& c1, 
					TstCalib2 const& c2) {
    std::array<double,3> aa;
    for(int channel=0; channel<3; channel++) {
      aa[channel] += c1.row(channel).dToE();
      if(c2.row(channel).status()!="OK") aa[channel] = -1.0;
    }
    return std::make_shared<DetData2>(aa);
  }

}
