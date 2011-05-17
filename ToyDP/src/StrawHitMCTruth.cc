// 
// This is a place to put additional information produced by HitMaker,
// like true drift distance, signal propagation time, etc.
//
// Original author Ivan Logashenko
//

// C++ includes
#include <ostream>

// Framework includes.
#include "cetlib/exception.h"

// Mu2e includes
#include "ToyDP/inc/StrawHitMCTruth.hh"

using namespace std;

namespace mu2e {

  // Print the information found in this hit.
  void StrawHitMCTruth::print( ostream& ost, bool doEndl ) const {

    ost << "Straw Hit MC Truth:"
        << " t0: "      << _t0
        << " drift time: "      << _driftTime
        << " drift distance: "      << _driftDistance
        << " distance to wire center: "     << _distanceToMid;

    if ( doEndl ){
      ost << endl;
    }
    
  }
  
} // namespace mu2e
