#include "DbExample/inc/DetData1Maker.hh"

namespace mu2e{

  DetData1::ptr_t DetData1Maker::fromStatic(double dtoe) {
    return std::make_shared<DetData1>(dtoe);
  }

  DetData1::ptr_t DetData1Maker::fromConfig(SimpleConfig const& sc) {
    return std::make_shared<DetData1>(sc.getDouble("DetData1.dtoe"));
  }

  DetData1::ptr_t DetData1Maker::fromDb(TstCalib1 const& t1) {
    float a = 0.0;
    for(auto const& r : t1.rows()) {
      a += r.dToE();
    }
    a /= t1.nrow();
    
    return std::make_shared<DetData1>(a);
  }

}
