//
// Extinction Monitor UCI Tof hit info plus possible additional information produced by HitMaker
//
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
