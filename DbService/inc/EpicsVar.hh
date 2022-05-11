#ifndef DbService_EpicsVar_hh
#define DbService_EpicsVar_hh

// This class holds one instance of one EPICS variable
// the value will be in of num_val, float_val or string_val

#include <ctime>
#include <string>
#include <vector>

namespace mu2e {

class EpicsVar {
 public:
  typedef std::vector<EpicsVar> EpicsVec;

  EpicsVar() {}
  EpicsVar(std::string const& csv);

  // the channel number internal to EPICS
  int8_t channel_id() const { return _channel_id; }
  // the time of the sample as a string, to us
  std::string const& stime() const { return _smpl_time_s; }
  // the time of the sample in epoch time
  std::time_t ttime() const { return _smpl_time_t; }
  // the nanosec past the intergral seconds
  int8_t nanosecs() const { return _nanosecs; }
  int8_t severity_id() const { return _severity_id; }
  int8_t status_id() const { return _status_id; }
  // PV value, if an int
  int num_val() const { return _num_val; }
  // PV value, if a float or double
  double float_val() const { return _float_val; }
  // PV value, if a string
  std::string const& str_val() const { return _str_val; }
  // this char says what the PV type is
  char datatype() const { return _datatype; }

  // the database row in csv
  std::string const& csv() const { return _csv; };

 private:
  std::string _csv;

  int8_t _channel_id;
  std::string _smpl_time_s;
  std::time_t _smpl_time_t;
  int8_t _nanosecs;
  int8_t _severity_id;
  int8_t _status_id;
  int _num_val;
  double _float_val;
  std::string _str_val;
  char _datatype;
};

}  // namespace mu2e

#endif
