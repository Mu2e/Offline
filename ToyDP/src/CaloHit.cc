//
// $Id: CaloHit.cc,v 1.2 2010/11/02 03:25:37 genser Exp $
// $Author: genser $
// $Date: 2010/11/02 03:25:37 $
//
// Original author Ivan Logashenko
//

// C++ includes
#include <ostream>

// Framework includes.
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "ToyDP/inc/CaloHit.hh"

using namespace std;

namespace mu2e {

  // Print the information found in this hit.
  void CaloHit::print( ostream& ost, bool doEndl ) const {

    ost << "Calorimeter Hit:   "
        << " RO id: "    << _roId
        << " time: "     << _time
        << " energyDep: "     << _energyDep;

    if ( doEndl ){
      ost << endl;
    }
    
  }
  
} // namespace mu2e
