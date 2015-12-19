//
// Original author G. Pezzullo
//

// C++ includes
#include <ostream>

// Framework includes.
#include "cetlib/exception.h"

// Mu2e includes
#include "RecoDataProducts/inc/CaloDigiPacked.hh"

using namespace std;

namespace mu2e {

  // Print the information found in this hit. TODO
  void CaloDigiPacked::print( ostream& ost, bool doEndl ) const {

    ost << " CaloDigiPacked output:   "
        << endl;


    if ( doEndl ){
      ost << endl;
    }

  }

} // namespace mu2e
