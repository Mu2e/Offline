#ifndef DbService_RunSubRun_hh
#define DbService_RunSubRun_hh

// holds info for one row in the subrun table

#include <string>
#include <vector>

namespace mu2e {

class RunSubRun {
 public:
  typedef std::vector<RunSubRun> RunSubRunVec;

  RunSubRun() :
      _run_number(0), _subrun_number(0), _n_events(0), _n_on_spill(0),
      _n_off_spill(0), _n_null(0), _min_ewt(0), _max_ewt(0),
      _start_time_unix(0), _stop_time_unix(0) {}

  int runNumber() const { return _run_number; }
  int subrunNumber() const { return _subrun_number; }
  long nEvents() const { return _n_events; }
  long nOnSpill() const { return _n_on_spill; }
  long nOffSpill() const { return _n_off_spill; }
  long nNull() const { return _n_null; }
  long minEwt() const { return _min_ewt; }
  long maxEwt() const { return _max_ewt; }
  int startTimeUnix() const { return _start_time_unix; }
  int stopTimeUnix() const { return _stop_time_unix; }
  const std::string& eventModeCounts() const { return _event_mode_counts; }
  const std::string& createdAt() const { return _created_at; }

  void setRunNumber(int run_number) { _run_number = run_number; }
  void setSubrunNumber(int subrun_number) { _subrun_number = subrun_number; }
  void setNEvents(long n_events) { _n_events = n_events; }
  void setNOnSpill(long n_on_spill) { _n_on_spill = n_on_spill; }
  void setNOffSpill(long n_off_spill) { _n_off_spill = n_off_spill; }
  void setNNull(long n_null) { _n_null = n_null; }
  void setMinEwt(long min_ewt) { _min_ewt = min_ewt; }
  void setMaxEwt(long max_ewt) { _max_ewt = max_ewt; }
  void setStartTimeUnix(int start_time_unix) {
    _start_time_unix = start_time_unix;
  }
  void setStopTimeUnix(int stop_time_unix) { _stop_time_unix = stop_time_unix; }
  void setEventModeCounts(const std::string& event_mode_counts) {
    _event_mode_counts = event_mode_counts;
  }
  void setCreatedAt(const std::string& created_at) { _created_at = created_at; }

 private:
  int _run_number;
  int _subrun_number;
  long _n_events;
  long _n_on_spill;
  long _n_off_spill;
  long _n_null;
  long _min_ewt;
  long _max_ewt;
  int _start_time_unix;
  int _stop_time_unix;
  std::string _event_mode_counts;  // jsonb stored as string
  std::string _created_at;
};

}  // namespace mu2e

#endif
