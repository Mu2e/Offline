//
// Object to allow exhaustively iterating over state permutations of a set of TrkStrawHits
//
// $Id: TrkStrawHitState.cc,v 1.1 2012/05/14 19:20:02 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/05/14 19:20:02 $
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

  TrkStrawHitState::TrkStrawHitState(TrkStrawHit* tsh) : _state(unknown), _tsh(tsh) {
    if(tsh->isActive()){
      if(tsh->ambig() < 0)
	_state = negambig;
      else if(tsh->ambig() < 0)
	_state = posambig;
      else
	_state = noambig;
    } else
      _state = inactive;
  }

  void
  TrkStrawHitState::setHitState(bool setactivity) const {
// the sign here is a result of the perverse convention in the BaBar
// code that a track displaced from the wire in the poAsitive 'U' direction has negative DOCA
    if(_tsh != 0){
      _tsh->setAmbigUpdate(true);
      switch (_state) {
	case inactive:
	  if(setactivity)_tsh->setActivity(false);
	  break;
	case unknown: default:
	  break;
	case noambig:
	  _tsh->setAmbig(0);
	  if(setactivity)_tsh->setActivity(true);
	  break;
	case negambig:
	  if(setactivity)_tsh->setActivity(true);
	  _tsh->setAmbig(1);
	  break;
	case posambig:
	  if(setactivity)_tsh->setActivity(true);
	  _tsh->setAmbig(-1);
	  break;
      }
    }
  }

  TrkStrawHitStateVector::TrkStrawHitStateVector() :  _state(0), _nstates(0) {
  }

  TrkStrawHitStateVector::TrkStrawHitStateVector(std::vector<TrkStrawHit*> const& tshv,TSHSSV const& allowed) : _allowed(allowed) {
    _tshsv.reserve(tshv.size());
    for(std::vector<TrkStrawHit*>::const_iterator itsh=tshv.begin();itsh != tshv.end(); ++itsh){
      _tshsv.push_back(TrkStrawHitState(*itsh));
    }
    _nstates = ipow(_allowed.size(),_tshsv.size());
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
  TrkStrawHitStateVector::setHitStates(bool setactivity) const {
    for(size_t itsh =0; itsh < _tshsv.size(); ++itsh){
	_tshsv[itsh].setHitState(setactivity);
    }
  }
  
  bool
  TrkStrawHitStateVector::increment() {
    if(_state < nStates()-1){
      ++_state;
      intToVect();
      return true;
    } else
      return false;
  }

  bool
  TrkStrawHitStateVector::decrement() {
    if(_state > 0){
      --_state;
      intToVect();
      return true;
    } else
      return false;
  }

  void
  TrkStrawHitStateVector::reset() {
    _state = 0;
    intToVect();
  }

  void
  TrkStrawHitStateVector::intToVect() {
// work from the highest bits to the lowest
// using a floating rep of the state
    unsigned state = _state;
    for(int itsh=_tshsv.size()-1;itsh>=0;--itsh){
      unsigned npow = ipow(_allowed.size(),(unsigned)itsh);
      unsigned istate = state/npow;
// check
      if(istate < _allowed.size()){
	_tshsv[itsh].setState(_allowed[istate]);
	state -= istate*npow;
      } else {
	throw cet::exception("RECO")<<"mu2e::TrkStrawHitStateVector: illegal state" << std::endl;
	state = 0.0;
      }
    }
// check for completeness
    if(state != 0)throw cet::exception("RECO")<<"mu2e::TrkStrawHitStateVector: illegal set of states" << std::endl;
  }

  void
  TrkStrawHitStateVector::vectToInt() {
    _state = 0;
    for(unsigned itsh=0;itsh<_tshsv.size();++itsh){
      TSHSSV::const_iterator ish = std::find(_allowed.begin(),_allowed.end(),_tshsv[itsh].state());
      if(ish == _allowed.end())throw cet::exception("RECO")<<"mu2e::TrkStrawHitStateVector: illegal state" << std::endl;
      _state += ((ish-_allowed.begin())/sizeof(TSHSSV::const_iterator))*ipow(_allowed.size(),itsh);
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

