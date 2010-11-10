//
// An enum-matched-to-names class for stopping codes from G4.
//
// $Id: StoppingCode.cc,v 1.1 2010/11/10 23:42:28 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/11/10 23:42:28 $
//
// Original author Rob Kutschke

#include <iostream>
#include <iomanip>

#include "ToyDP/inc/StoppingCode.hh"

#include <boost/static_assert.hpp>

using namespace std;

namespace mu2e {
  
  const char* StoppingCode::_name[] = { STOPPINGCODE_NAMES };

  BOOST_STATIC_ASSERT(sizeof(StoppingCode::_name)/sizeof(char*) == StoppingCode::lastEnum);

  void StoppingCode::printAll( std::ostream& ost){
    ost << "List of stopping codes from G4: " << endl;
    for ( int i=0; i<lastEnum; ++i){
      ost << setw(2) << i << " " << _name[i] << std::endl;
    }
  }

}
