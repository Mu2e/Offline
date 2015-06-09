//
// Object to allow exhaustively iterating over state permutations of a set of TrkStrawHits
//
// $Id: HitState.cc,v 1.2 2012/05/22 21:35:42 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/05/22 21:35:42 $
//
// Original author David Brown, LBNL
//
#include "KalmanTests/inc/PanelAmbigStructs.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "cetlib/coded_exception.h"
#include <iostream>
#include <algorithm>

namespace mu2e 
{
  namespace PanelAmbig {

    HitState::HitState(int ambig, bool active) {
      if(!active)
	_state = inactive;
      else {
	if(ambig > 0)
	  _state = posambig;
	else if(ambig < 0)
	  _state = negambig;
	else
	  _state = noambig;
      }
    }

    HitState::HitState(const TrkStrawHit* tsh)  {
      if(tsh->isActive()){
	if(tsh->ambig() < 0)
	  _state = negambig;
	else if(tsh->ambig() > 0)
	  _state = posambig;
	else
	  _state = noambig;
      } else
	_state = inactive;
    }

    void HitState::setHitState(TrkStrawHit* tsh) const {
      if(tsh != 0){
	switch (_state) {
	  default:
	    break;
	  case noambig:
	    tsh->setActivity(true);
	    tsh->setAmbig(0);
	    break;
	  case negambig:
	    tsh->setActivity(true);
	    tsh->setAmbig(-1);
	    break;
	  case posambig:
	    tsh->setActivity(true);
	    tsh->setAmbig(1);
	    break;
	  case inactive:
	    tsh->setActivity(false);
	    break;
	}
      }
    }

    void PanelState::setHitStates(TSHV& tshv) const {
      for(size_t ihit=0;ihit<tshv.size(); ++ihit)
	hitState(ihit).setHitState(tshv[ihit]);
    }

    PanelResult::PanelResult(PanelState const& state) : _state(state),
    _chisq(-1.0), _status(-1), _hpat(null), _statechange(0) {
// classify the state.  If any hits are opposite, it's opposite.  Otherwise if at least 2
// are on the same side, they are same.  Otherwise it is null
      for(size_t ihit=0;ihit<_state._nhits;++ihit){
	for(size_t jhit=ihit+1;jhit<_state._nhits;++jhit){
	  int hprod = static_cast<int>(_state.hitState(ihit)._state) * 
	    static_cast<int>(_state.hitState(jhit)._state);
	    if(hprod < 0){
	      _hpat = opposite;
	      break;
	    } else if(hprod > 0)
	      _hpat = same;
	}
	if(_hpat == opposite)break;// double break
      }
    }

    TSHUInfo::TSHUInfo(const TrkStrawHit* tsh,CLHEP::Hep3Vector const& udir, HepPoint const& uorigin) : _use(free), _hstate(tsh) {
      // find wire position at POCA	
      HepPoint wpos = tsh->hitTraj()->position(tsh->hitLen());
      CLHEP::Hep3Vector wposv(wpos.x(),wpos.y(),wpos.z());
      // translate WRT origin and project
      CLHEP::Hep3Vector dstraw = wpos-uorigin;
      _upos = udir.dot(dstraw);
      _wcpos = udir.dot(wposv);
      // note that the t0 component of the error scales coherently between the hits, so here we use only the intrinsic error
      _uerr = tsh->hitErr();
      _uwt = 1.0/(_uerr*_uerr);
      _dr = tsh->driftRadius();
      _dv = tsh->driftVelocity();
      _ambig = tsh->ambig();
      _active = tsh->isActive();
      _index = tsh->index();
    }

  } // PanelAmbig namespace
} // mu2e namespace

