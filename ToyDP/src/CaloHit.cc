
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

    ost << "Calorimeter Hit:"
        << " RO id: "    << _roId
        << " time: "     << _time
        << " eDep: "     << _energyDep;

    if ( doEndl ){
      ost << endl;
    }
    
  }
  
} // namespace mu2e
