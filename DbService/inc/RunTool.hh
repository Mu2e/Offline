#ifndef DbService_RunTool_hh
#define DbService_RunTool_hh

// This class takes a user request, queries the run_info offline
// replica database and returns information about the run configuration

#include "Offline/DbService/inc/DbReader.hh"
#include "Offline/DbService/inc/RunInfo.hh"
#include "Offline/DbService/inc/RunSelect.hh"

#include <string>

namespace mu2e {

class RunTool {
 public:
  enum class FlagType { run = 0, transition = 1 };

  RunTool();

  // fill argument with the names of all flag types
  std::map<int, std::string> flags(FlagType ftype);
  // Both flag types (run, transition) as a single JSON object -- string
  // return, same SWIG-compatibility reason as runJson(). Lives on the class
  // rather than in runTool_main.cc so it's reachable from the Python/SWIG
  // binding too, not just the CLI.
  std::string flagsJson();
  RunInfo::RunVec listRuns(const RunSelect& runsel, bool configs = false,
                           bool transitions = false, bool subruns = false);
  void printRun(const RunInfo& info);
  // Same content as printRun, as a JSON object (serialized text, not a
  // nlohmann::json return type -- this class is SWIG-wrapped and
  // nlohmann::json isn't a SWIG-recognized type, same reason
  // RunInfo::dbTables2/3() and RunConfig::dbTables2/3() return std::string).
  std::string runJson(const RunInfo& info);

 private:
  // read the database via query_engine
  DbReader _reader;
};

}  // namespace mu2e

#endif
