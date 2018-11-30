#ifndef DbExample_DetData1_hh
#define DbExample_DetData1_hh

// time dependant detector data 1
// with dependence one db table
// this code is in an object directory, last in link order

#include <string>
#include <memory>
#include "DbExample/inc/ConditionsEntity2.hh"

namespace mu2e {

  class DetData1 : public ConditionsEntity2 {
  public:
    typedef std::shared_ptr<DetData1> ptr_t;
    DetData1():DetData1(0) {}
    DetData1(double data):_name("DetData1"),_data(data) {}
    double getData() const { return _data;}
    std::string const& name() const { return _name; }
  private:
    std::string _name;
    double _data;
  };

};

#endif
