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

    HitState::HitState(TrkStrawHit* tsh) : _state(ignore){
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
	  case ignore : default:
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
	}
      }
    }

    PanelState::PanelState() : _nhits(0)  {
      for(size_t ihit=0;ihit<MAXNHITS; ++ihit)
	_hstate[ihit] = HitState(HitState::ignore);
    }

    void PanelState::setHitStates(TSHV& tshv) const {
      for(size_t ihit=0;ihit<tshv.size(); ++ihit)
	_hstate[ihit].setHitState(tshv[ihit]);
    }

  } // PanelAmbig namespace
} // mu2e namespace

