//
// Crystal hit info plus possible additional information produced by HitMaker
//
// $Id: CaloCrystalOnlyHit.cc,v 1.2 2011/06/07 21:32:21 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/07 21:32:21 $
//

// Mu2e includes
#include "MCDataProducts/inc/CaloCrystalOnlyHit.hh"

using namespace std;

namespace mu2e {

  void CaloCrystalOnlyHit::setEnergyDep(double energy) {
    _energyDep = energy;
    return;
  }

  // Print the information found in this hit.
  void CaloCrystalOnlyHit::print( ostream& ost, bool doEndl ) const {

    ost << "Calorimeter Crystal Hit MC:"
        << " crystal id: "  << _crystalId
        << " time "         << _time
        << " energyDep: "   << _energyDep;

    if ( doEndl ){
      ost << endl;
    }

  }

} // namespace mu2e
