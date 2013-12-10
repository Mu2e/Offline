//
// StrawHitletSequence is a time-ordered sequence of StrawHitlets
//
// $Id: StrawHitletSequence.cc,v 1.2 2013/12/10 01:32:51 brownd Exp $
// $Author: brownd $
// $Date: 2013/12/10 01:32:51 $
//
// Original author David Brown, LBNL
//
// mu2e includes
#include "TrackerMC/inc/StrawHitletSequence.hh"
#include "cetlib/exception.h"

using namespace std;

namespace mu2e {
  StrawHitletSequence::StrawHitletSequence() : _strawIndex(0), _end(StrawEnd::unknown)
  {}

  StrawHitletSequence::StrawHitletSequence(StrawHitlet const& hitlet) :
    _strawIndex(hitlet.strawIndex()), _end(hitlet.strawEnd())
  {
    insert(hitlet);
  }

  StrawHitletSequence::StrawHitletSequence(StrawIndex const& index, StrawEnd end) :
    _strawIndex(index), _end(end)
  {}
  
  StrawHitletSequence::StrawHitletSequence(StrawHitletSequence const& other) : 
    _strawIndex(other._strawIndex), _end(other._end), _hlist(other._hlist) {}

  StrawHitletSequence& StrawHitletSequence::operator =(StrawHitletSequence const& other) {
    if(&other != this){
      _strawIndex = other._strawIndex;
      _end = other._end;
      _hlist = other._hlist;
    }
    return *this;
  }
  // insert a new hitlet.  This is the only non-trivial function
  HitletList::iterator StrawHitletSequence::insert(StrawHitlet const& hitlet) {
    HitletList::iterator retval = _hlist.end();
    if(hitlet.type() == StrawHitlet::unknown){
      throw cet::exception("SIM") 
	<< "mu2e::StrawHitletSequence: tried to add unknown hitlet type" 
	<< endl;
      return retval;
    }
    // make sure the straw and end are the same
    if(!_hlist.empty() && (hitlet.strawIndex() != strawIndex() 
      || hitlet.strawEnd() != strawEnd())){
      throw cet::exception("SIM") 
	<< "mu2e::StrawHitletSequence: tried to add hitlet from a different straw/end to a sequence" 
	<< endl;
      return retval;
    }
    if(_hlist.empty()){
      _hlist.push_front(hitlet);
      _strawIndex = hitlet.strawIndex();
      _end = hitlet.strawEnd();
      retval = _hlist.begin();
    } else {
      // loop over the contents and insert in the correct place
      HitletList::iterator ibefore = _hlist.begin();
      while(ibefore != _hlist.end() && ibefore->time() < hitlet.time())
	++ibefore;
      retval = _hlist.insert(ibefore,hitlet);
    }
    return retval;
  }
}

