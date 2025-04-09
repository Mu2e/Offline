//
// Object to allow exhaustively iterating over state permutations of a set of TrkStrawHits
//
//
// Original author David Brown, LBNL
//
#include "Offline/TrkReco/inc/PanelAmbigStructs.hh"
#include "cetlib_except/coded_exception.h"
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

    bool HitState::setHitState(TrkStrawHit* tsh) const {
      bool oldactive = tsh->isActive();
      int oldambig = tsh->ambig();
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
      return !(oldactive == tsh->isActive() && oldambig == tsh->ambig());
    }

    PanelResult::PanelResult(PanelState const& state) : _state(state),
    _chisq(-1.0), _status(-1), _hpat(null), _statechange(0) {
      // classify the state.  If any hits are opposite, it's opposite.  Otherwise if at least 2
      // are on the same side, they are same.  Otherwise it is null
      for(size_t ihit=0;ihit<_state.size();++ihit){
        for(size_t jhit=ihit+1;jhit<_state.size();++jhit){
          int hprod = static_cast<int>(_state[ihit]._state) *
            static_cast<int>(_state[jhit]._state);
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

    unsigned ipow(unsigned base, unsigned exp) {
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

    bool setHitStates(PanelState const& pstate, TrkStrawHitVector& hits) {
      bool retval(false); // assume nothing changes
      if(pstate.size() != hits.size())
        throw cet::exception("RECO")<<"mu2e::PanelAmbigResolver: state size doesn't match" << std::endl;
      for(size_t ihit=0;ihit<hits.size(); ++ihit)
        retval |= pstate[ihit].setHitState(hits[ihit]);
      return retval;
    }

  } // PanelAmbig namespace
} // mu2e namespace

