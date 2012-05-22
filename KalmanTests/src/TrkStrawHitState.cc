//
// Object to allow exhaustively iterating over state permutations of a set of TrkStrawHits
//
// $Id: TrkStrawHitState.cc,v 1.2 2012/05/22 21:35:42 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/05/22 21:35:42 $
//
// Original author David Brown, LBNL
//
#include "KalmanTests/inc/TrkStrawHitState.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "cetlib/coded_exception.h"
#include <iostream>
#include <algorithm>

namespace mu2e 
{

  TrkStrawHitState::TrkStrawHitState(TrkStrawHit* tsh) : _state(ignore), _tsh(tsh) {
    if(tsh->isActive()){
      if(tsh->ambig() > 0)
	_state = negambig;
      else if(tsh->ambig() < 0)
	_state = posambig;
      else
	_state = noambig;
    } else
      _state = inactive;
  }

  void
  TrkStrawHitState::setHitState() const {
// the sign here is a result of the perverse convention in the BaBar
// code that a track displaced from the wire in the poAsitive 'U' direction has negative DOCA
    if(_tsh != 0){
      _tsh->setAmbigUpdate(true);
      switch (_state) {
	case inactive:
	  _tsh->setActivity(false);
	  break;
	case ignore : default:
	  break;
	case noambig:
	  _tsh->setActivity(true);
	  _tsh->setAmbig(0);
	  break;
	case negambig:
	  _tsh->setActivity(true);
	  _tsh->setAmbig(1);
	  break;
	case posambig:
	  _tsh->setActivity(true);
	  _tsh->setAmbig(-1);
	  break;
      }
    }
  }

  TrkStrawHitStateVector::TrkStrawHitStateVector() :  _nstates(0) {
  }

  TrkStrawHitStateVector::TrkStrawHitStateVector(TSHV const& tshv,TSHSSV const& allowed) : _allowed(allowed) {
    _tshsv.reserve(tshv.size());
    for(TSHV::const_iterator itsh=tshv.begin();itsh != tshv.end(); ++itsh){
// only record those hits whose states are currently allowed.  The others will be held frozen
      TrkStrawHitState tshs(*itsh);
      TSHSSV::const_iterator ish = std::find(_allowed.begin(),_allowed.end(),tshs.state());
      if(ish != _allowed.end())
	_tshsv.push_back(tshs);
      else
	_tshsv.push_back(TrkStrawHitState(TrkStrawHitState::ignore,*itsh));
    }
    vectToInt();
  }

  TrkStrawHitStateVector::TrkStrawHitStateVector(TrkStrawHitStateVector const& other) : _tshsv(other._tshsv), _allowed(other._allowed), _state(other._state),
    _nstates(other._nstates)
  {}
  
  TrkStrawHitStateVector&
  TrkStrawHitStateVector::operator =(TrkStrawHitStateVector const& other) {
    if(&other != this){
      _tshsv = other._tshsv;
      _allowed = other._allowed;
      _nstates = other._nstates;
      _state = other._state;
    }
    return *this;
  }

  void
  TrkStrawHitStateVector::setHitStates() const {
    for(size_t itsh =0; itsh < _tshsv.size(); ++itsh){
	_tshsv[itsh].setHitState();
    }
  }
  
  bool
  TrkStrawHitStateVector::increment() {
    if(incrementState()){
      intToVect();
      return true;
    } else
      return false;
  }

  bool
  TrkStrawHitStateVector::decrement() {
    if(decrementState()){
      intToVect();
      return true;
    } else
      return false;
  }

  void
  TrkStrawHitStateVector::reset() {
    resetState();
    intToVect();
  }

  bool
  TrkStrawHitStateVector::incrementState() {
    bool retval(false);
// increment the first possible state
    for(size_t itsh=0;itsh<_tshsv.size();++itsh){
      if(_tshsv[itsh].state() != TrkStrawHitState::ignore){
	if(_state[itsh] < (int)(_allowed.size()-1)){
	  ++_state[itsh];
	  retval = true;
	  break;
	} else {
	  _state[itsh] = 0;
	}
      }
    }
    return retval;
  }

  bool
  TrkStrawHitStateVector::decrementState() {
    bool retval(false);
// decrement the first possible state
    for(size_t itsh=0;itsh<_tshsv.size();++itsh){
      if(_tshsv[itsh].state() != TrkStrawHitState::ignore){
	if(_state[itsh] > 0){
	  --_state[itsh];
	  retval = true;
	  break;
	} else {
	  _state[itsh] = _allowed.size()-1;
	}
      }
    }
    return retval;
  }

  void 
  TrkStrawHitStateVector::resetState() {
    for(size_t itsh=0;itsh<_tshsv.size();++itsh){
      if(_tshsv[itsh].state() != TrkStrawHitState::ignore){
	_state[itsh] = 0;
      }
    }
  }

  void
  TrkStrawHitStateVector::intToVect() {
    for(size_t itsh=0;itsh<_tshsv.size();++itsh){
      if(_tshsv[itsh].state() != TrkStrawHitState::ignore){
	_tshsv[itsh].setState(_allowed[_state[itsh]]);
      }
    }
  }

  void
  TrkStrawHitStateVector::vectToInt() {
    _nstates = 1;
    for(size_t itsh=0;itsh<_tshsv.size();++itsh){
      if(_tshsv[itsh].state() != TrkStrawHitState::ignore){
	TSHSSV::const_iterator ish = std::find(_allowed.begin(),_allowed.end(),_tshsv[itsh].state());
	if(ish == _allowed.end())throw cet::exception("RECO")<<"mu2e::TrkStrawHitStateVector: illegal state" << std::endl;
	size_t state = (ish-_allowed.begin());
	_state.push_back(state);
	_nstates *= _allowed.size();
      } else {
	_state.push_back(-1);
      }
    }
  }

  unsigned
  TrkStrawHitStateVector::ipow(unsigned base, unsigned exp)
  {
    unsigned result = 1;
    while (exp)
    {
      if (exp & 1)
	result *= base;
      exp >>= 1;
      base *= base;
    }
    return result;
  }
}

