#ifndef DbService_EpicsTool_hh
#define DbService_EpicsTool_hh

// This class takes a user request, queries the EPICS offline replica database
// and returns info about epics variables

#include "Offline/DbService/inc/DbReader.hh"
#include "Offline/DbService/inc/EpicsVar.hh"
#include <ctime>
#include <string>
#include <vector>

namespace mu2e {

class EpicsTool {
 public:
  typedef std::vector<std::string> StringVec;

  EpicsTool();

  // fill argument with the names of all channels
  int names(StringVec& names);
  // return a vector of EPICS samples of channel "name"
  // and consistent with time. Time can be a single time
  // 2022-05-05T12:45:02-05:00
  // and the cloest sample before the time is returned.
  // a range of times
  // 2022-05-05T12:45:02-05:00/2022-05-05T13:40:02-05:00
  // where samples between the times is returned
  // or a float (1.5 for example), and the result will be sample
  // between this many days ago and now
  EpicsVar::EpicsVec get(std::string const& name, std::string const& time);

 private:
  // turn time string into a epoch time
  std::time_t parseTimeTZ(std::string const& stime) const;
  // read the database via query_engine
  DbReader _reader;
};

}  // namespace mu2e

#endif
