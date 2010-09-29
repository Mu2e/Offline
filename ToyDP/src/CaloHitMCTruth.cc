// 
// This is a place to put additional information produced by HitMaker,
// like true drift distance, signal propagation time, etc.
//
// Original author Ivan Logashenko
//

// C++ includes
#include <ostream>

// Framework includes.
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "ToyDP/inc/CaloHitMCTruth.hh"

using namespace std;

namespace mu2e {

  // Print the information found in this hit.
  void CaloHitMCTruth::print( ostream& ost, bool doEndl ) const {

    ost << "Calorimeter Hit MC Truth:"
        << " crystal id: "      << _roId
        << " energy dep: "      << _energyDep;

    if ( doEndl ){
      ost << endl;
    }
    
  }
  
} // namespace mu2e
