// 
// This is a place to put additional information produced by HitMaker,
// like true drift distance, signal propagation time, etc.
//
// $Id: CaloCrystalHitMCTruth.cc,v 1.1 2010/10/28 20:43:58 genser Exp $
// $Author: genser $
// $Date: 2010/10/28 20:43:58 $
//

// C++ includes
#include <ostream>

// Framework includes.
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "ToyDP/inc/CaloCrystalHitMCTruth.hh"

using namespace std;

namespace mu2e {

  // Print the information found in this hit.
  void CaloCrystalHitMCTruth::print( ostream& ost, bool doEndl ) const {

    ost << "Calorimeter Crystal Hit MC Truth:"
        << " crystal id: "  << _crystalId
        << " time "         << _time
        << " energy dep: "  << _energyDep
        << " used roids: "  << _numberOfROIdsUsed
        << " crystal roids:";
    for (size_t i=0; i!=_roIds.size(); ++i) {
      ost  << " " << i << ": " <<_roIds[i];
    }

    if ( doEndl ){
      ost << endl;
    }
    
  }
  
} // namespace mu2e
