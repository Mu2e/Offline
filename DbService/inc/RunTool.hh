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
  enum class FlagType { run = 0, transition = 1, cause = 2 };

  RunTool();

  // fill argument with the names of all channels
  std::map<int, std::string> flags(FlagType ftype);
  RunInfo::RunVec listRuns(const RunSelect& runsel, bool conditions = false,
                           bool transitions = false);
  void printRun(const RunInfo& info, bool longp = false);

 private:
  // read the database via query_engine
  DbReader _reader;
};

}  // namespace mu2e

#endif
