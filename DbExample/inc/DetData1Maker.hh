#ifndef DbExample_DetData1Maker_hh
#define DbExample_DetData1Maker_hh

// make time dependant - detector data 1
// with dependence on one db table
// this file and cc will be in the detector code area

#include "ConfigTools/inc/SimpleConfig.hh"
#include "DbExample/inc/DetData1.hh"
#include "DbTables/inc/TstCalib1.hh"

namespace mu2e {

  class DetData1Maker {
  public:
    DetData1::ptr_t fromStatic(double x);
    DetData1::ptr_t fromConfig(SimpleConfig const& sc);
    DetData1::ptr_t fromDb(TstCalib1 const& t1);
  };

};

#endif
