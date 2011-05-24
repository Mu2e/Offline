//
// Status information about running G4 for one event.
//
// $Id: StatusG4.cc,v 1.1 2011/05/24 17:16:44 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/24 17:16:44 $
//
// Original author Rob Kutschke
//

#include <iostream>

// Mu2e includes
#include "MCDataProducts/inc/StatusG4.hh"

using namespace std;

namespace mu2e{

  void StatusG4::swap( StatusG4& rhs ){
    std::swap( _status,               rhs._status);
    std::swap( _nG4Tracks,            rhs._nG4Tracks);
    std::swap( _overflowSimParticles, rhs._overflowSimParticles);
    std::swap( _nKilledStepLimit,     rhs._nKilledStepLimit);
    std::swap( _cpuTime,              rhs._cpuTime);
    std::swap( _realTime,             rhs._realTime);
  }

  void StatusG4::print ( ostream& ost, bool newLine ) const {
    ost << "G4 status: " << _status
        << "; Number G4 Tracks: "
        << _nG4Tracks << ";";
    if ( _overflowSimParticles ){
      ost << " SimParticleCollection overflowed;";
    }
    ost << "   Number killed by too many steps: " << _nKilledStepLimit
        << ";";
    if ( newLine ){
      ost << endl;
    }
  }
}
