#ifndef DbService_RunInfo_hh
#define DbService_RunInfo_hh

// holds info for one run for runTool

#include <string>
#include <vector>

namespace mu2e {

class RunInfo {
 public:
  typedef std::vector<RunInfo> RunVec;

  RunInfo(std::string csv, std::string conditions, std::string transitions);

  int runNumber() const { return _run_number; }
  int runType() const { return _run_type; }
  int conditionId() const { return _condition_id; }
  int artdaqPartition() const { return _artdaq_partition; }
  const std::string& hostName() const { return _host_name; }
  const std::string& configurationName() const { return _configuration_name; }
  int configurationVersion() const { return _configuration_version; }
  const std::string& contextName() const { return _context_name; }
  int contextVersion() const { return _context_version; }
  const std::string& triggerTableName() const { return _trigger_table_name; }
  int triggerTableVersion() const { return _trigger_table_version; }
  const std::string& onlineSoftwareVersion() const {
    return _online_software_version;
  }
  const std::string& commitTime() const { return _commit_time; }
  const std::string& shifterNote() const { return _shifter_note; }
  const std::string& conditions() const { return _conditions; }
  const std::string& transitions() const { return _transitions; }

  const std::string& csv() const { return _csv; }

 private:
  std::string _csv;

  int _run_number;
  int _run_type;
  int _condition_id;
  int _artdaq_partition;
  std::string _host_name;
  std::string _configuration_name;
  int _configuration_version;
  std::string _context_name;
  int _context_version;
  std::string _trigger_table_name;
  int _trigger_table_version;
  std::string _online_software_version;
  std::string _commit_time;
  std::string _shifter_note;

  std::string _conditions;
  std::string _transitions;
};

}  // namespace mu2e

#endif
