//
// First version of a hit as described by Mu2e-doc-900.
//
// $Id: StrawHit.cc,v 1.1 2011/05/24 17:16:43 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/24 17:16:43 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <ostream>

// Framework includes.
#include "cetlib_except/exception.h"

// Mu2e includes
#include "RecoDataProducts/inc/StrawHit.hh"

using namespace std;

namespace mu2e {

  // Print the information found in this hit.
  void StrawHit::print( ostream& ost, bool doEndl ) const {

    ost << "traw Hit:"
        << " idx: "      << _strawIndex
        << " time: "     << _time
        << " dt: "       << _dt
        << " eDep: "     << _energyDep;

    if ( doEndl ){
      ost << endl;
    }

  }

} // namespace mu2e
