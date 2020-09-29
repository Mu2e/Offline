// C++ includes
#include <iostream>

// Mu2e includes
#include "TrackerConditions/inc/TrackerStatus.hh"
#include <sstream>

using namespace std;

namespace mu2e {

  void TrackerStatus::print( ostream& out) const{
    for( auto const& estat :  _estatus) {
	  out << "Tracker Element with Id " << estat.sid_
	    << " level " << estat.mask_.levelName()  
	    << " status " << estat.status_ << std::endl;
    }
  }
  
} // namespace mu2e
