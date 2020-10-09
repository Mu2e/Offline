//
// First version of a hit as described by Mu2e-doc-900.
//
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

    ost << " StrawHit:"
        << " idx: "      << _strawId
        << " cal time "     << _time[StrawEnd::cal]
        << " HV time "     << _time[StrawEnd::hv]
        << " cal TOT "     << _tot[StrawEnd::cal]
        << " HV TOT "     << _tot[StrawEnd::hv]
        << " eDep: "     << _energyDep;

    if ( doEndl ){
      ost << endl;
    }

  }

} // namespace mu2e
