// 
// Crystal hit info plus possible additional information produced by HitMaker
//
// $Id: CaloCrystalOnlyHit.cc,v 1.1 2010/12/13 06:10:33 logash Exp $
// $Author: logash $
// $Date: 2010/12/13 06:10:33 $
//

// C++ includes
#include <ostream>

// Framework includes.
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "ToyDP/inc/CaloCrystalOnlyHit.hh"
#include "ToyDP/inc/DPIndex.hh"

using namespace std;

namespace mu2e {

  void CaloCrystalOnlyHit::setEnergyDep(double energy) {
    _energyDep = energy;
    return;
  }

  // Print the information found in this hit.
  void CaloCrystalOnlyHit::print( ostream& ost, bool doEndl ) const {

    ost << "Calorimeter Crystal Hit MC:"
        << " crystal id: "  << _crystalId
        << " time "         << _time
        << " energyDep: "   << _energyDep;
 
    if ( doEndl ){
      ost << endl;
    }
    
  }
  
} // namespace mu2e
