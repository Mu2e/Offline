//
// $Id: CaloHit.cc,v 1.4 2011/05/18 02:27:19 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:19 $
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
