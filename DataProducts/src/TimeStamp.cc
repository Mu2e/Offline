// A presistent time stamp class that can represent times as seconds from the start of the
// unix epoch, in UTC, to the 32-bit unsigned epoch rollover in 2106.
//
// Original author Rob Kutschke

#include "Offline/DataProducts/inc/TimeStamp.hh"
#include <iomanip>
#include <ostream>

namespace mu2e {

  // Print time, in UTC, as a formatted string, including the time zone specifier.
  std::ostream& operator<<(std::ostream& os,
			   TimeStamp const& ts ){
    auto t = ts.get();
    os << std::put_time(gmtime(&t),"%FT%T%z");
    return os;
  }

}
