//
// Original author Stefano Roberto Soleti
//

// C++ includes
#include <ostream>

// Framework includes.
#include "cetlib/exception.h"

// Mu2e includes
#include "RecoDataProducts/inc/CaloDigi.hh"

using namespace std;

namespace mu2e {

  // Print the information found in this hit. TODO
  void CaloDigi::print( ostream& ost, bool doEndl ) const {

    ost << " CaloDigi output:   "
        << endl;


    if ( doEndl ){
      ost << endl;
    }

  }

} // namespace mu2e
