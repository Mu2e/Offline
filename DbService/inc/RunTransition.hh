#ifndef DbService_RunTransition_hh
#define DbService_RunTransition_hh

// holds info for one row in the run_transition table

#include <string>
#include <vector>

namespace mu2e {

class RunTransition {
 public:
  typedef std::vector<RunTransition> RunTransitionVec;

  RunTransition() : _run_number(0), _type_id(0), _cause_id(0) {}

  int runNumber() const { return _run_number; }
  int typeId() const { return _type_id; }
  int causeId() const { return _cause_id; }
  const std::string& transitionTime() const { return _transition_time; }

  void setRunNumber(int run_number) { _run_number = run_number; }
  void setTypeId(int type_id) { _type_id = type_id; }
  void setCauseId(int cause_id) { _cause_id = cause_id; }
  void setTransitionTime(const std::string& transition_time) {
    _transition_time = transition_time;
  }

 private:
  int _run_number;
  int _type_id;
  int _cause_id;
  std::string _transition_time;
};

}  // namespace mu2e

#endif
