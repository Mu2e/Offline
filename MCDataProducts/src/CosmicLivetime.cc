//
// Original author Stefano Roberto Soleti
//

// C++ includes
#include <ostream>

// Framework includes.
#include "cetlib_except/exception.h"

// Mu2e includes
#include "MCDataProducts/inc/CosmicLivetime.hh"

using namespace std;

namespace mu2e {

  // Print the information found in this hit.
  void CosmicLivetime::print( ostream& ost, bool doEndl ) const {

    ost << "Cosmic livetime: " << _livetime << " s" << endl
        << " primaries: " << _primaries << endl
        << " area: "      << _area << " m^2" << endl
        << " lowE: "       << _lowE << " GeV" << endl
        << " highE: "  << _highE << " GeV" << endl
        << " fluxConstant: "     << _fluxConstant;

    if ( doEndl ){
      ost << endl;
    }

  }

} // namespace mu2e
