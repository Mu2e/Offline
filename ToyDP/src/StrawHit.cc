// 
// First version of a hit as described by Mu2e-doc-900.
//
// $Id: StrawHit.cc,v 1.2 2010/08/18 23:14:03 logash Exp $
// $Author: logash $
// $Date: 2010/08/18 23:14:03 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <ostream>

// Framework includes.
#include "FWCore/Utilities/interface/Exception.h"

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
