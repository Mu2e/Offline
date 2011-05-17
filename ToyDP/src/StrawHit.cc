// 
// First version of a hit as described by Mu2e-doc-900.
//
// $Id: StrawHit.cc,v 1.3 2011/05/17 15:36:01 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:36:01 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <ostream>

// Framework includes.
#include "cetlib/exception.h"

// Mu2e includes
#include "ToyDP/inc/StrawHit.hh"

using namespace std;

namespace mu2e {

  // Print the information found in this hit.
  void StrawHit::print( ostream& ost, bool doEndl ) const {

    ost << "traw Hit:"
        << " idx: "      << _strawIndex
        << " time: "     << _time
        << " dt: "       << _dt
        << " eDep: "     << _energyDep;

    if ( doEndl ){
      ost << endl;
    }
    
  }
  
} // namespace mu2e
