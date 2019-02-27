#ifndef DbExample_DetData3_hh
#define DbExample_DetData3_hh

// time dependant detector data 3
// with dependence on only two conditions objects
// not the db tables directly
// this code is in an object directory, last in link order

#include <string>
#include <array>
#include <memory>
#include "DbExample/inc/ConditionsEntity2.hh"

namespace mu2e {

  class DetData3 : public ConditionsEntity2 {
  public:
    typedef std::shared_ptr<DetData3> ptr_t;
    DetData3():DetData3(std::array<double,3>()) {}
    DetData3(std::array<double,3> data):_name("DetData3"),_data(data) {}
    std::array<double,3> getData() const { return _data;}
    std::string const& name() const { return _name; }
  private:
    std::string _name;
    std::array<double,3> _data;
  };

};

#endif
