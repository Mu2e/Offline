#ifndef DbExample_DetData3Maker_hh
#define DbExample_DetData3Maker_hh

// make time dependant - detector data 3 
// with dependence on only two conditions objects
// not the db tables directly
// this file and cc will be in the detector code area


#include <vector>
#include "ConfigTools/inc/SimpleConfig.hh"
#include "DbExample/inc/DetData1.hh"
#include "DbExample/inc/DetData2.hh"
#include "DbExample/inc/DetData3.hh"


namespace mu2e {

  class DetData3Maker {
  public:
    DetData3::ptr_t fromConfig(SimpleConfig const& sc);
    DetData3::ptr_t fromStatic(std::vector<double> const& vv);
    DetData3::ptr_t fromDb(DetData1 const& d1, DetData2 const& d2);
    
  };

};

#endif
