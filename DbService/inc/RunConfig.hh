#ifndef DbService_RunConfig_hh
#define DbService_RunConfig_hh

// holds info for one row in the run_info config table

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
  std::string dbTables3(bool qjson = false) const;

  // Walk the entire JSON tree in _settings, collect all (cid, name) pairs
  // from every "DBServiceCIDTable" array found at any depth.
  // If json=true  -> returns a JSON list of {"cid":..., "name":...} objects.
  // If json=false -> returns one "cid name" pair per line.
  std::string dbTables2(bool qjson = false) const;

  // Walk the entire JSON tree in _settings for the "TriggerConfigTable" node
  // and summarize the trigger menu it carries: the table's identity, and for
  // each trigger path its bit, enabled flag, and the prescale and output
  // streams of every event mode.
  // If qjson=true  -> returns a JSON object with "found":true, in the fixed
  //                   schema triggerMenuNotFound() documents.
  // If qjson=false -> returns a human-readable block: an identity line, then
  //                   the paths sorted by bit with one prescale column per
  //                   event mode, then the logger streams.
  // Returns an empty string if this config carries no TriggerConfigTable
  // (same convention as dbTables2/3 above).
  std::string triggerMenu(bool qjson = false) const;

  // The object triggerMenu(true) returns, with "found":false and every other
  // field at its empty default -- for a run with no TRG config at all, or one
  // whose TRG config holds no TriggerConfigTable. Lives here, beside the
  // builder, so the schema's key list is written in exactly one place;
  // RunInfo::triggerMenu() is the caller.
  static std::string triggerMenuNotFound(int run_number);

 private:
  int _run_number;
  std::string _subsystem;
  std::string _settings;  // jsonb stored as string
  std::string _create_time;
  int _version;
};

}  // namespace mu2e

#endif
