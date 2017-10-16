#ifndef Alignment_IoV_HH
#define Alignment_IoV_HH
// IoV.hh
// This is the header file for the IoV class, instances of which are
// Intervals of Validity.
// This is a placeholder IoV scheme until one is decided on by
// the collaboration.
// This implementation consists of a start time and end time,
// represented as doubles.  That may change later.

#include <limits>

namespace mu2e {

  class IoV {
  public:
    IoV( double start, double end ) :
    _startTime(start),
      _endTime(end) {}

    double start() const {return _startTime;}
    double end()   const {return _endTime;}

    double _maxtime = std::numeric_limits<double>::max() - 10.0;  // allow for something...?

    static bool isValid(double aTime) { return !(aTime < 0 || aTime > (std::numeric_limits<double>::max() -10.0) ); }
    bool contains(double aTime) const { return (aTime >= _startTime && aTime <= _endTime); }

  private:

    double _startTime;
    double _endTime;
  };

} // end namespace mu2e

#endif //  Alignment_IoV_HH
