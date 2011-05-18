//
// First version of a hit as described by Mu2e-doc-900.
//
// $Id: StrawHit.cc,v 1.4 2011/05/18 02:27:19 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:19 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <ostream>

// Framework includes.
#include "cetlib/exception.h"

// Mu2e includes
#include "ToyDP/inc/StrawHit.hh"

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
