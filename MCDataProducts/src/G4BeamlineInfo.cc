// C++ includes
#include <ostream>

// Framework includes.
#include "cetlib_except/exception.h"

// Mu2e includes
#include "MCDataProducts/inc/G4BeamlineInfo.hh"

using namespace std;

namespace mu2e {

  // Print the information found in this hit.
  void G4BeamlineInfo::print( ostream& ost, bool doEndl ) const {

    ost << "G4Beamline extra data:"
        << " eventId: "      << _event_id
        << " trackId: "      << _track_id
        << " weight: "       << _weight
        << " time: "         << _time;

    if ( doEndl ){
      ost << endl;
    }

  }

  void G4BeamlineInfo::swap(G4BeamlineInfo & a) {
    G4BeamlineInfo tmp(a);
    a = *this;
    *this = tmp;
  }

} // namespace mu2e
