#ifndef DbService_EpicsVar_hh
#define DbService_EpicsVar_hh

// This class holds one instance of one EPICS variable
// the value will be in of num_val, float_val or string_val

#include <ctime>
#include <string>
#include <vector>
#include <variant>

namespace mu2e {

class EpicsVar {
 public:
  typedef std::vector<EpicsVar> EpicsVec;
  typedef std::variant<int,double,std::string> variVar;

  EpicsVar(std::string const& csv, int64_t channel_id,
           std::string const& smpl_time_s, std::time_t smpl_time_t,
           int64_t nanosecs, int64_t severity_id,
           int64_t status_id, variVar const& value) :
    _csv(csv),_channel_id(channel_id),_smpl_time_s(smpl_time_s),
    _smpl_time_t(smpl_time_t),_nanosecs(nanosecs),_severity_id(severity_id),
    _status_id(status_id),_value(value) {}

  // the channel number internal to EPICS
  int64_t channel_id() const { return _channel_id; }
  // the time of the sample as a string, to us
  std::string const& stime() const { return _smpl_time_s; }
  // the time of the sample in epoch time
  std::time_t ttime() const { return _smpl_time_t; }
  // the nanosec past the integral seconds
  int64_t nanosecs() const { return _nanosecs; }
  int64_t severity_id() const { return _severity_id; }
  int64_t status_id() const { return _status_id; }
  // PV value held, with its type, as std::variant
  variVar value() const { return _value; }
  // for convenience, testing and querying variant
  bool isInt() const { return std::holds_alternative<int>(_value); }
  int getInt() const { return std::get<int>(_value); }
  bool isDouble() const { return std::holds_alternative<double>(_value); }
  double getDouble() const { return std::get<double>(_value); }
  bool isString() const { return std::holds_alternative<std::string>(_value); }
  std::string getString() const { return std::get<std::string>(_value); }

  // the database row in csv
  std::string const& csv() const { return _csv; };

 private:
  std::string _csv;

  int64_t _channel_id;
  std::string _smpl_time_s;
  std::time_t _smpl_time_t;
  int64_t _nanosecs;
  int64_t _severity_id;
  int64_t _status_id;
  variVar _value;
};

}  // namespace mu2e

#endif
