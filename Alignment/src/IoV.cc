// IoV.cc
// Definitions file for the IoV (Interval of Validity) class
// 
// David Norvil Brown, UofL, Oct 2017

#include "Alignment/inc/IoV.hh"

namespace mu2e {

  IoV::IoV( const IoV& rhs ) {
    _startTime = rhs.start();
    _endTime = rhs.end();
  } // end of copy c-tor

  std::ostream& operator<<(std::ostream& os, const IoV& rhs ) {
    os << "IoV: (" << rhs.start() << " => " << rhs.end() << ") \n";
    return os;
  } // end of output operator


} // end namespace mu2e
