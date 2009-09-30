//
// A persistable G4 hit representing the intersection of a simulated
// particle with straw gas.
//
//
// $Id: StrawMCHit.cc,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke

// Mu2e incldues
#include "ToyDP/inc/StrawMCHit.hh"

using namespace std;

namespace mu2e {


  // Should switch these two to the standard idiom for 
  // copy c'tor and assignment.
  StrawMCHit::StrawMCHit(const StrawMCHit& right):
    _trackId(right._trackId),
    _strawIndex(right._strawIndex),
    _edep(right._edep),
    _time(right._time),
    _position(right._position),
    _momentum(right._momentum){
  }

  const StrawMCHit& StrawMCHit::operator=(const StrawMCHit& right){
    _trackId    = right._trackId;
    _strawIndex = right._strawIndex;
    _edep       = right._edep;
    _time       = right._time;
    _position   = right._position;
    _momentum   = right._momentum;
    return *this;
  }

  void StrawMCHit::print( ostream& ost ) const {

    ost << "  trackId: "        << _trackId 
	<< "  strawIndex: "     << _strawIndex
	<< "  energy deposit: " << _edep
	<< "  position: "       << _position
	<< "  momentum: "       << _momentum
	<< "  time: "           << _time
	<< endl;
  }

} // namespace mu2e
