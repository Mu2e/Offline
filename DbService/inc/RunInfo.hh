#ifndef DbService_RunInfo_hh
#define DbService_RunInfo_hh

// holds info for one run, including related records from multiple tables

#include "Offline/DbService/inc/RunConfig.hh"
#include "Offline/DbService/inc/RunRecord.hh"
#include "Offline/DbService/inc/RunSubRun.hh"
#include "Offline/DbService/inc/RunTransition.hh"

#include <vector>

namespace mu2e {

class RunInfo {
 public:
  typedef std::vector<RunInfo> RunVec;

  RunInfo() {}

  // Access to the main run record
  const RunRecord& runRecord() const { return _run_record; }
  RunRecord& runRecord() { return _run_record; }
  void setRunRecord(const RunRecord& run_record) { _run_record = run_record; }

  // Access to configuration records for this run
  const RunConfig::RunConfigVec& configs() const { return _configs; }
  RunConfig::RunConfigVec& configs() { return _configs; }
  void setConfigs(const RunConfig::RunConfigVec& configs) {
    _configs = configs;
  }
  void addConfig(const RunConfig& config) { _configs.push_back(config); }

  // Access to transition records for this run
  const RunTransition::RunTransitionVec& transitions() const {
    return _transitions;
  }
  RunTransition::RunTransitionVec& transitions() { return _transitions; }
  void setTransitions(const RunTransition::RunTransitionVec& transitions) {
    _transitions = transitions;
  }
  void addTransition(const RunTransition& transition) {
    _transitions.push_back(transition);
  }

  // Access to subrun records for this run
  const RunSubRun::RunSubRunVec& subruns() const { return _subruns; }
  RunSubRun::RunSubRunVec& subruns() { return _subruns; }
  void setSubruns(const RunSubRun::RunSubRunVec& subruns) {
    _subruns = subruns;
  }
  void addSubrun(const RunSubRun& subrun) { _subruns.push_back(subrun); }

  // Convenience accessor for run number
  int runNumber() const { return _run_record.runNumber(); }

  // Check if run is done and return status code
  // Returns:
  //   0 = still running
  //   1 = user stop (transitions 0,1,6,7)
  //   2 = error stop (transition 2)
  //   3 = probably crash end (timeout exceeded)
  // Throws exception if subruns or transitions are empty
  int isDone(int timeout_seconds = 3600) const;

 private:
  RunRecord _run_record;
  RunConfig::RunConfigVec _configs;
  RunTransition::RunTransitionVec _transitions;
  RunSubRun::RunSubRunVec _subruns;
};

}  // namespace mu2e

#endif
