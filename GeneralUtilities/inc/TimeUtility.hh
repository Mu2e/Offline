#ifndef GeneralUtilities_TimeUtility_hh
#define GeneralUtilities_TimeUtility_hh

// a time conversion method needed a few places

#include <ctime>
#include <string>

namespace mu2e {

class TimeUtility {
public:
  // allowed formats:
  //  2018-10-12 (first second of this date, UTC)
  //  2018-10-12T08:58:26 (UTC)
  // in the future, the above should probably default to FNAL local time
  //  2018-10-12T08:58:26-05:00
  //  2018-10-12T08:58:26.792518-05:00
  // and the above without the T
  //  2018-10-12 08:58:26.792518-05:00
  // and the above with an alternative timezone format
  //  2018-10-12T08:58:26.792518-0500
  // will return time in ttime, return non-zero if parse error
  static int parseTimeTZ(std::string const& stime, std::time_t& ttime);
  //  2018-10-12 08:58:26.792518-05:00 to 2018-10-12T08:58:26-05:00
  static std::string reformat1(std::string const& stime);
};
} // namespace mu2e
#endif
