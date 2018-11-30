#ifndef DbExample_DetData2_hh
#define DbExample_DetData2_hh

// time dependant detector data 2
// with dependence on two db tables
// this code is in an object directory, last in link order

#include <string>
#include <array>
#include <memory>
#include "DbExample/inc/ConditionsEntity2.hh"

namespace mu2e {

  class DetData2 : public ConditionsEntity2 {
  public:
    typedef std::shared_ptr<DetData2> ptr_t;
    DetData2():DetData2(std::array<double,3>()) {}
    DetData2(std::array<double,3> data):_name("DetData2"),_data(data) {}
    std::array<double,3> getData() const { return _data;}
    std::string const& name() const { return _name; }
  private:
    std::string _name;
    std::array<double,3> _data;
  };

};

#endif
