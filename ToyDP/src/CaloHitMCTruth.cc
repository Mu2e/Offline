// 
// This is a place to put additional information produced by HitMaker,
//
// $Id: CaloHitMCTruth.cc,v 1.3 2010/11/11 21:17:49 genser Exp $
// $Author: genser $
// $Date: 2010/11/11 21:17:49 $
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

    ost << "Calorimeter Hit MC:"
        << " RO id: "      << _roId
        << " time: "       << _time
        << " energyDep: "  << _energyDep
        << " charge: "     << _charged;

    if ( doEndl ){
      ost << endl;
    }
    
  }
  
} // namespace mu2e
