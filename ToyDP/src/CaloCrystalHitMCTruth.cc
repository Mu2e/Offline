// 
// Crystal hit info plus possible additional information produced by HitMaker
//
// $Id: CaloCrystalHitMCTruth.cc,v 1.3 2010/11/11 21:16:53 genser Exp $
// $Author: genser $
// $Date: 2010/11/11 21:16:53 $
//

// C++ includes
#include <ostream>

// Framework includes.
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "ToyDP/inc/CaloCrystalHitMCTruth.hh"
#include "ToyDP/inc/CaloHitMCTruth.hh"
#include "ToyDP/inc/DPIndex.hh"

using namespace std;

namespace mu2e {

  // Print the information found in this hit.
  void CaloCrystalHitMCTruth::print( ostream& ost, bool doEndl ) const {

    ost << "Calorimeter Crystal Hit MC:"
        << " crystal id: "  << _crystalId
        << " time "         << _time
        << " energyDep: "  << _energyDep;
 
    if ( doEndl ){
      ost << endl;
    }
    
  }
  
} // namespace mu2e
