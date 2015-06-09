//
// Object to allow exhaustively iterating over state permutations of a set of TrkStrawHits
//
// $Id: HitState.hh,v 1.2 2012/05/22 21:35:42 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/05/22 21:35:42 $
//
#ifndef mu2e_PanelAmbig_PanelAmbigStructs_HH
#define mu2e_PanelAmbig_PanelAmbigStructs_HH
#include "KalmanTests/inc/threevec.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "Rtypes.h"
#include <cstddef>
#include <vector>
#include <cassert>

class HepPoint;

namespace mu2e {
  class TrkStrawHit;
  namespace PanelAmbig {
// struct to describe the ambiguity and activity state of a single TrkStrawHit
    struct HitState {
      enum TSHState {negambig=-1,noambig,posambig,inactive};
      HitState(TSHState state=negambig) : _state(state){}
      HitState(int ambig, bool active);
      // set state from an existing strawhit
      HitState(const TrkStrawHit* tsh);
      // set TrkStrawHit to this state
      void setHitState(TrkStrawHit* tsh) const;
      // comparators
      bool operator == (HitState const& other) const { return _state == other._state; }
      bool operator != (HitState const& other) const { return !operator==(other); }
      bool operator < (HitState const& other ) const { return _state < other._state; }
      TSHState _state; //
    };
    typedef std::vector<TrkStrawHit*> TSHV;
    typedef std::vector<HitState> HSV;

    // class for root to hold Hit States per panel.  Genreflex can't
    // handle unions, arrays, ..., so these are hard coded FIXME!!! 
  #define MAXNHITS 6
    struct PanelHitState {
      HitState _hit0, _hit1, _hit2, _hit3, _hit4, _hit5;
      HitState& hitState(size_t ish) {
	static HitState def(HitState::inactive);
	assert(ish <= MAXNHITS);
	switch(ish) {
	  case 0:
	    return _hit0;
	  case 1:
	    return _hit1;
	  case 2:
	    return _hit2;
	  case 3:
	    return _hit3;
	  case 4:
	    return _hit4;
	  case 5:
	    return _hit5;
	  default:
	    return def;
	}
	return def;
      }
      HitState const& hitState(size_t ish) const {
	PanelHitState* ncthis = const_cast<PanelHitState*>(this);
	return ncthis->hitState(ish);
      }
    };

    // compact description of the state of an entire panel
    struct PanelState {
      PanelHitState _pstate;
      unsigned _nhits;
      PanelState() : _nhits(0) {}
      // construct from a vector of states
      PanelState(HSV const& hsv) : _nhits(hsv.size()) {
	assert(hsv.size() <=MAXNHITS);
	for(size_t ish=0;ish < hsv.size(); ++ish)
	  _pstate.hitState(ish) = hsv[ish];
      }
      // construct from a vector of hits.  This assumes all hits are from the same panel!!!
      PanelState(TSHV const& tshv) : _nhits(tshv.size()) {
	assert(tshv.size() <=MAXNHITS);
	for(size_t ish=0;ish < tshv.size(); ++ish)
	  _pstate.hitState(ish) = HitState(tshv[ish]);
      }
      // set the state of the referenced TrkStrawHits.
      void setHitStates(TSHV& hits) const;
      // accessors
      HitState& hitState(size_t ish) { return _pstate.hitState(ish); }
      HitState const& hitState(size_t ish) const { return _pstate.hitState(ish); }
    };

    typedef std::vector<PanelState> PSV;

   // class to store the 1-d projection of a hit.
    // this direction (u) is perpendicular to the wire and to the track as it hits the panel
    struct TSHUInfo {
      enum Use {free=0,fixed,unused}; // free hits are allowed change.  Fixed are not, but are still used in chisquared.  unused hits are ignored
      // construct from a TrkStrawHit, the U direction and the U origin in mu2e coordinates
      TSHUInfo(const TrkStrawHit* tsh,CLHEP::Hep3Vector const& udir, HepPoint const& uorigin);
      TSHUInfo() : _use(free), _upos(0.0), _wcpos(0.0), _uerr(-1.0), _uwt(0.0), _dr(0.0), _dv(0.0), _ambig(0), _active(-1) {}
      Use _use; // use status of this hit
      HitState _hstate; // initial hitstate
      Float_t _upos; // hit position in u relative to the projected track position
      Float_t _wcpos; // wire center position in detector coordinate
      Float_t _uerr; // hit error projected in u direction 
      Float_t _uwt; // cache of the hit weight = 1/error^2
      Float_t _dr; // drift radius
      Float_t _dv; // drift velocity
      Int_t _ambig; // TrkStrawHit nominal ambiguity
      Int_t _active; // TrkStrawHit nominal activity
      UInt_t _index; // index into TSH array
    };
    typedef std::vector<TSHUInfo> TSHUIV;

    // MC truth hit U projection information, for diagnostics
    struct TSHMCUInfo {
      Float_t _dr; // true drift radius
      Int_t _ambig;
      TSHMCUInfo() : _dr(-1000.0),_ambig(-10) {}
    };
    typedef std::vector<TSHMCUInfo> TSHMCUIV;

    // struct to hold hit and other information for a single panel.
    // It also holds the result of the optimization
    struct PanelInfo {
      PanelInfo() : _tupos(0.0), _tuerr(0.0), _tuwt(0.0), _nused(0), _nfree(0) {}
      // information from the track fit, excluding the hits from panel, projected into this panel's measurement direction
      threevec _udir;  // u direction in detector coordinates
      Float_t _tupos; // hit position projected perpendicular to track direction WRT wire
      Float_t _tuerr; // track error projected into U coordinate
      Float_t _tuwt; // cache of track weight = 1/error^2
      UInt_t _nused; // # of hits used in chisquared
      UInt_t _nfree; // # of hits whose state is allowed to vary
      // information from all the hits
      TSHUIV _uinfo;
    };

    // struct to hold the result of optimizing the panel parameters.
    struct PanelResult {
    // classify the hit patterns in terms of having hits on the same or opposite sides
    // null is when
      enum HitPattern{null=0,same,opposite};
      PanelResult(PanelState const& state);
      PanelResult() : _chisq(-1.0), _status(-1), _hpat(null) , _statechange(0) {}
      PanelState _state; // ambiguity/activity state of each hit in this panel
      Float_t _chisq; // chisquared of this state, including all penalties
      CLHEP::HepVector _delta; // change in parameters at optimum; U is parameter 0, deltaT is parameter 1
      CLHEP::HepSymMatrix _dcov; // covariance of parameter changes
      Int_t _status; // status of matrix inversion
      HitPattern _hpat;
      UInt_t _statechange; // bit pattern of hit state changes
    };

    typedef std::vector<PanelResult> PRV;
  }
}
#endif
