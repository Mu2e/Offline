//
// $Id: CaloHit.cc,v 1.3 2011/05/17 15:36:01 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:36:01 $
//
// Original author Ivan Logashenko
//

// C++ includes
#include <ostream>

// Framework includes.
#include "cetlib/exception.h"

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
