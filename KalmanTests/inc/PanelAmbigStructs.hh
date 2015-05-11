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
// class to describe the ambiguity and activity state of a single TrkStrawHit
    class HitState {
      public:
	enum TSHState {negambig=-1,noambig,posambig,inactive,ignore,nstates};
	// set state explicitly.  
	HitState(TSHState state=ignore) : _state(state){}
	// set state from an existing strawhit
	HitState(TrkStrawHit* tsh);
	// return the state as an integer
	TSHState state() const { return _state; }
	// total number of states
	static unsigned nStates() { return nstates; }
	void setState(TSHState state) { _state = state; }
	// set TrkStrawHit tot his state
	void setHitState(TrkStrawHit* tsh) const;
	// comparators
	bool operator == (HitState const& other) const { return _state == other._state; }
	bool operator != (HitState const& other) const { return !operator==(other); }
	bool operator < (HitState const& other ) const { return _state < other._state; }
      private:
	TSHState _state;
    };
    typedef std::vector<TrkStrawHit*> TSHV;
    typedef std::vector<HitState> HSV;

    // compact dscription of the state of an entire panel
    struct PanelState {
#define MAXNHITS 4
      HitState _hstate[MAXNHITS];
      size_t _nhits;
      PanelState();
      // construct from a vector of states
      PanelState(HSV const& hsv) : _nhits(hsv.size()) {
	assert(hsv.size() <=MAXNHITS);
	for(size_t ish=0;ish < hsv.size(); ++ish)
	  _hstate[ish] = hsv[ish];
      }
      // construct from a vector of hits.  This assumes all hits are from the same panel!!!
      PanelState(TSHV const& tshv) : _nhits(tshv.size()) {
	assert(tshv.size() <=MAXNHITS);
	for(size_t ish=0;ish < tshv.size(); ++ish)
	  _hstate[ish] = HitState(tshv[ish]);
      }
      // set the state of the referenced TrkStrawHits.
      void setHitStates(TSHV& hits) const;

      // accessors
      HitState const& state(size_t ish) const {
	return _hstate[4];
      }
    };

    typedef std::vector<PanelState> PSV;

   // class to store the 1-d projection of a hit.
    // this direction (u) is perpendicular to the wire and to the track as it hits the panel
    struct TSHUInfo {
      // construct from a TrkStrawHit, the U direction and the U origin in mu2e coordinates
      TSHUInfo(const TrkStrawHit* tsh,CLHEP::Hep3Vector const& udir, HepPoint const& uorigin);
      TSHUInfo() : _tsh(0), _upos(0.0), _wcpos(0.0), _uerr(-1.0), _uwt(0.0), _dr(0.0), _dv(0.0), _ambig(0), _active(-1) {}
      const TrkStrawHit* _tsh; // reference the straw hit, to allow consistency checking
      Float_t _upos; // hit position in u relative to the projected track position
      Float_t _wcpos; // wire center WRT track 
      Float_t _uerr; // hit error projected in u direction 
      Float_t _uwt; // cache of the hit weight = 1/error^2
      Float_t _dr; // drift radius
      Float_t _dv; // drift velocity
      Int_t _ambig; // TrkStrawHit nominal ambiguity
      Int_t _active; // TrkStrawHit nominal activity
    };
    typedef std::vector<TSHUInfo> TSHUIV;

    // MC truth hit U projection information, for diagnostics
    struct TSHMCUInfo {
      Float_t _upos;
      Int_t _ambig;
    };
    typedef std::vector<TSHMCUInfo> TSHMCUIV;

    // struct to hold hit and other information for a single panel.
    // It also holds the result of the optimization
    struct PanelInfo {
      PanelInfo() : _tupos(0.0), _tuerr(0.0), _tuwt(0.0) {}
      // information from the track fit, excluding the hits from panel, projected into this panel's measurement direction
      threevec _udir;  // u direction in detector coordinates
      Float_t _tupos; // hit position projected perpendicular to track direction WRT wire
      Float_t _tuerr; // track error projected into U coordinate
      Float_t _tuwt; // cache of track weight = 1/error^2
      // information from all the hits
      TSHUIV _uinfo;
    };

    // struct to hold the result of optimizing the panel parameters.
    struct PanelResult {
      PanelResult(PanelState const& state) : _state(state),
      _chisq(-1.0), _status(-1) {}
      PanelResult() : _chisq(-1.0), _status(-1) {}
      PanelState _state; // ambiguity/activity state of each hit in this panel
      Float_t _chisq; // chisquared of this state, including all penalties
      CLHEP::HepVector _delta; // change in parameters at optimum; U is parameter 0, deltaT is parameter 1
      CLHEP::HepSymMatrix _dcov; // covariance of parameter changes
      Int_t _status; // status of matrix inversion
    };

    typedef std::vector<PanelResult> PRV;
  }
}
#endif
