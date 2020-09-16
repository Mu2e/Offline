//
// Status information about running G4 for one event.
//
//
// Original author Rob Kutschke
//

#include "MCDataProducts/inc/StatusG4.hh"

#include <iostream>

using namespace std;

namespace mu2e{

  void StatusG4::swap( StatusG4& rhs ){
    std::swap( _status,               rhs._status);
    std::swap( _nG4Tracks,            rhs._nG4Tracks);
    std::swap( _overflowSimParticles, rhs._overflowSimParticles);
    std::swap( _nKilledStepLimit,     rhs._nKilledStepLimit);
    std::swap( _nKilledByFieldPropagator, rhs._nKilledByFieldPropagator);
    std::swap( _cpuTime,              rhs._cpuTime);
    std::swap( _realTime,             rhs._realTime);
  }

  void StatusG4::print ( ostream& ost, bool newLine ) const {

    if ( _status == undefined ){
      ost << "G4 status: undefined";
    } else{


      ost << "G4 status: " << _status
          << "; Number G4 Tracks: "
          << _nG4Tracks << ";";
      if ( _overflowSimParticles ){
        ost << " SimParticleCollection overflowed;";
      }
      ost << "   Number killed by too many steps: " << _nKilledStepLimit
          << "   Number killed by Field Propagator: " << _nKilledByFieldPropagator
          << ";";
    }

    if ( newLine ){
      ost << endl;
    }

  } // end StatusG4::print

  // Add the contents of another StatusG4 object to this one.
  // Used by event mixing to form a summary from multiple events.
  void StatusG4::add( StatusG4 const& rhs ){

    // Set status to the worst of the two contributing statuses.
    if ( _status == undefined ){
      _status = rhs._status;
    } else{
      _status = ( _status > rhs._status ) ? _status : rhs._status;
    }

    _nG4Tracks            += rhs._nG4Tracks;
    _overflowSimParticles += rhs._overflowSimParticles;
    _nKilledStepLimit     += rhs._nKilledStepLimit;
    _nKilledByFieldPropagator += rhs._nKilledByFieldPropagator;
    _cpuTime              += rhs._cpuTime;
    _realTime             += rhs._realTime;

  }

} // end namespace mu2e

