#ifndef DbService_EpicsTool_hh
#define DbService_EpicsTool_hh

// This class takes a user request, queries the EPICS offline replica database
// and returns epics variables names and samples content

#include "Offline/DbService/inc/DbReader.hh"
#include "Offline/DbService/inc/EpicsVar.hh"
#include "Offline/GeneralUtilities/inc/StringVec.hh"
#include <ctime>
#include <string>
#include <vector>

namespace mu2e {

class EpicsTool {
 public:

  EpicsTool();

  // fill argument with the names of all channels
  int names(StringVec& names);
  // return a vector of EPICS samples of channel "name"
  // and consistent with time or daysAgo. Time can be a single time
  // 2022-05-05T12:45:02-05:00
  // and the closest sample before the time is returned.
  // or a range of times
  // 2022-05-05T12:45:02-05:00/2022-05-05T13:40:02-05:00
  // where samples between the times is returned.
  // If time is empty and daysAgo is non-zero, then
  // the result will be samples between this many days ago and now
  EpicsVar::EpicsVec get(std::string const& name, std::string const& time,
                         float daysAgo = 0.0);

 private:
  // read the database via query_engine
  DbReader _reader;
};

}  // namespace mu2e

#endif
