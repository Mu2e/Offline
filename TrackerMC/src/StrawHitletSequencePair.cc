//
// StrawHitletSequencePair
//
// $Id: StrawHitletSequencePair.cc,v 1.1 2013/12/07 19:51:42 brownd Exp $
// $Author: brownd $
// $Date: 2013/12/07 19:51:42 $
//
// Original author David Brown, LBNL
//
// mu2e includes
#include "TrackerMC/inc/StrawHitletSequencePair.hh"
#include "cetlib/exception.h"

namespace mu2e {  
  StrawHitletSequencePair::StrawHitletSequencePair() {}
  StrawHitletSequencePair::StrawHitletSequencePair(StrawIndex index) : 
    _plus(index,StrawEnd::plus),
    _minus(index,StrawEnd::minus)
  {}

  StrawHitletSequencePair::StrawHitletSequencePair(StrawHitletSequencePair const& other) :
  _plus(other._plus), _minus(other._minus)
  {}

  StrawHitletSequencePair& StrawHitletSequencePair::operator =(StrawHitletSequencePair const& other) {
    if(&other != this){
      _plus = other._plus;
      _minus = other._minus;
    }
    return *this;
  }

  void StrawHitletSequencePair::insert(StrawHitletPair const& hpair) {
    if(hpair[0].strawEnd() != StrawEnd::minus || 
	hpair[1].strawEnd() != StrawEnd::plus ||
	hpair[0].strawIndex() != hpair[1].strawIndex())	

      throw cet::exception("SIM") 
	<< "mu2e::StrawHitletSequence: tried to add inconsistent hitlet pair to sequence";

// maybe could use move symanatics here?  FIXME!!
    _minus = hpair[0];
    _plus = hpair[1];

  }
  
}

