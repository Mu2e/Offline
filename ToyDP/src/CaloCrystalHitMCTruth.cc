// 
// Crystal hit info plus possible additional information produced by HitMaker
//
// $Id: CaloCrystalHitMCTruth.cc,v 1.4 2010/11/12 21:43:54 genser Exp $
// $Author: genser $
// $Date: 2010/11/12 21:43:54 $
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

  void CaloCrystalHitMCTruth::setEnergyDep(double energy) {
    _energyDep = energy;
    return;
  }

  // Print the information found in this hit.
  void CaloCrystalHitMCTruth::print( ostream& ost, bool doEndl ) const {

    ost << "Calorimeter Crystal Hit MC:"
        << " crystal id: "  << _crystalId
        << " time "         << _time
        << " energyDep: "   << _energyDep;
 
    if ( doEndl ){
      ost << endl;
    }
    
  }
  
} // namespace mu2e
