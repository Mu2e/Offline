#ifndef DbExample_DetData2Maker_hh
#define DbExample_DetData2Maker_hh

// make time dependant - detector data 2
// with dependence on one db table
// this file and cc will be in the detector code area


#include <vector>
#include "ConfigTools/inc/SimpleConfig.hh"
#include "DbExample/inc/DetData2.hh"
#include "DbTables/inc/TstCalib1.hh"
#include "DbTables/inc/TstCalib2.hh"

namespace mu2e {

  class DetData2Maker {
  public:
    DetData2::ptr_t fromStatic(std::vector<double> const& vv);
    DetData2::ptr_t fromConfig(SimpleConfig const& sc);
    DetData2::ptr_t fromDb(TstCalib1 const& c1, TstCalib2 const& c2);
  };

};

#endif
