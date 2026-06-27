#ifndef DbService_RunConfig_hh
#define DbService_RunConfig_hh

// holds info for one row in the config table

#include <string>
#include <vector>

namespace mu2e {

class RunConfig {
 public:
  typedef std::vector<RunConfig> RunConfigVec;

  RunConfig() : _run_number(0), _version(0) {}

  int runNumber() const { return _run_number; }
  const std::string& subsystem() const { return _subsystem; }
  const std::string& settings() const { return _settings; }
  const std::string& createTime() const { return _create_time; }
  int version() const { return _version; }

  void setRunNumber(int run_number) { _run_number = run_number; }
  void setSubsystem(const std::string& subsystem) { _subsystem = subsystem; }
  void setSettings(const std::string& settings) { _settings = settings; }
  void setCreateTime(const std::string& create_time) {
    _create_time = create_time;
  }
  void setVersion(int version) { _version = version; }

  // Walk the entire JSON tree in _settings, collect all key-value pairs
  // from every "DBServiceTables" dictionary found at any depth.
  // If json=true  -> returns a JSON object string with all merged pairs.
  // If json=false -> returns one VALUE per line (flat list).
  std::string dbTables3(bool json = false) const;

 private:
  int _run_number;
  std::string _subsystem;
  std::string _settings;  // jsonb stored as string
  std::string _create_time;
  int _version;
};

}  // namespace mu2e

#endif
