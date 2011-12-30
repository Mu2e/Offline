//
// Extinction Monitor UCI Tof hit info plus possible additional information produced by HitMaker
//
// $Id: ExtMonUCITofHit.cc,v 1.1 2011/12/30 20:31:46 youzy Exp $
// $Author: youzy $
// $Date: 2011/12/30 20:31:46 $
//

// Mu2e includes
#include "RecoDataProducts/inc/ExtMonUCITofHit.hh"

using namespace std;

namespace mu2e {

  void ExtMonUCITofHit::setEnergyDep(double energy) {
    _energyDep = energy;
    return;
  }

  // Print the information found in this hit.
  void ExtMonUCITofHit::print( ostream& ost, bool doEndl ) const {

    ost << "ExtMonUCI Hit :"
        << " station id: "  << _stationId
        << " segment id: "  << _segmentId
        << " time "         << _time
        << " energyDep: "   << _energyDep;

    if ( doEndl ){
      ost << endl;
    }

  }

} // namespace mu2e
