#ifndef Alignment_IoV_HH
#define Alignment_IoV_HH
// IoV.hh
// This is the header file for the IoV class, instances of which are
// Intervals of Validity.
// This is a placeholder IoV scheme until one is decided on by
// the collaboration.
// This implementation consists of a start time and end time,
// represented as unsigned longs.  That may change later.

#include <limits>
#include <ostream>

namespace mu2e {

  class IoV {
  public:
    IoV( const unsigned long& start, const unsigned long& end ) :
    _startTime(start),
      _endTime(end) {}

    IoV() :
      _startTime(0.0),
      _endTime( std::numeric_limits<unsigned long>::max() - 10 ) {}

    // copy c-tor
    IoV( const IoV& rhs );


    double start() const {return _startTime;}
    double end()   const {return _endTime;}

    static constexpr unsigned long maxtime = 
      std::numeric_limits<unsigned long>::max() - 10;// allow for something...?

    static bool isValid(const unsigned long& aTime) 
      { return !(aTime < 0 || aTime > (std::numeric_limits<unsigned long>::max() - 10) ); }

    bool contains(const unsigned long&  aTime) const 
      { return (aTime >= _startTime && aTime <= _endTime); }

    void updateEndTime( const unsigned long& aTime ) { _endTime = aTime; }

    IoV operator+ ( const IoV& rhs );

  private:

    unsigned long _startTime;
    unsigned long _endTime;

  }; // end of class def

  std::ostream& operator<<(std::ostream& os, const IoV& rhs);

} // end namespace mu2e

#endif //  Alignment_IoV_HH
