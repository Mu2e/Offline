#ifndef DbService_RunRecord_hh
#define DbService_RunRecord_hh

// holds info for one row in the runs table

#include <string>
#include <vector>

namespace mu2e {

class RunRecord {
 public:
  typedef std::vector<RunRecord> RunRecordVec;

  RunRecord() : _run_number(0), _run_type_id(0) {}

  int runNumber() const { return _run_number; }
  const std::string& comment() const { return _comment; }
  const std::string& createTime() const { return _create_time; }
  int runTypeId() const { return _run_type_id; }

  void setRunNumber(int run_number) { _run_number = run_number; }
  void setComment(const std::string& comment) { _comment = comment; }
  void setCreateTime(const std::string& create_time) {
    _create_time = create_time;
  }
  void setRunTypeId(int run_type_id) { _run_type_id = run_type_id; }

 private:
  int _run_number;
  std::string _comment;
  std::string _create_time;
  int _run_type_id;
};

}  // namespace mu2e

#endif
