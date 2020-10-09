//
//
// Original author Ivan Logashenko
//

// C++ includes
#include <ostream>

// Framework includes.
#include "cetlib_except/exception.h"

// Mu2e includes
#include "RecoDataProducts/inc/CaloHit.hh"

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
