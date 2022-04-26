#ifndef DataProducts_TimeStamp_hh
#define DataProducts_TimeStamp_hh
//
// A presistent time stamp class that represents a time, as seconds from the start of the
// unix epoch.  Mu2e plans to store times in UTC but this is not enforced by this class.
// The range of representable times is from the start of the unix epoch to the 32-bit 
// unsigned epoch rollover on Feb 27, 2106.
//
// Notes:
//  1) On all of our machines
//       - the C/C++ type time_t is int64_t.
//       - These wall clock functions return time as seconds since the unix epoch in UTC
//            - c++: std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
//            - c:   time(0);
//       - Since they are in UTC, there is no issue with DST transitions.
//       - with 8 bytes, epoch rollover is many times the lifetime of the universe away
//  2) If we truncate to int64_t, the epoch rollover is in 2038.
//  3) If we truncate to uint64_t, the epoch follover is in 2106.  We choose this.
//
// Original author Rob Kutschke
//

#include <cstdint>
#include <ctime>
#include <iosfwd>

namespace mu2e {

  class  TimeStamp {

  public:

    typedef uint32_t TimeStamp_t;

    TimeStamp(){}

    TimeStamp( time_t time ):
      time_(static_cast<TimeStamp_t>(time)){
    }

    time_t get()        const { return static_cast<time_t>( time_); }
    time_t operator()() const { return get(); }

  private:

     TimeStamp_t time_ = 0;

  };

  // Print time as a formatted string in UTC, including the time zone specifier.
  // For other formats or time zone options, call std::put_time directly.
  // Fixme: do we want the default to be Fermilab time?  This might be supported in C++20?
  std::ostream& operator<<(std::ostream& os,
			   TimeStamp const& ts );

}
#endif
