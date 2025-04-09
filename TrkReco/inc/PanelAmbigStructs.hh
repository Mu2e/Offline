//
// Object to allow exhaustively iterating over state permutations of a set of TrkStrawHits
//
// Original author: David Brown (LBNL), 2012
//
//
#ifndef mu2e_PanelAmbig_PanelAmbigStructs_HH
#define mu2e_PanelAmbig_PanelAmbigStructs_HH
#include "Offline/DataProducts/inc/GenVector.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "Rtypes.h"
#include "Offline/BTrkData/inc/TrkStrawHit.hh"
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
      bool setHitState(TrkStrawHit* tsh) const;
      // comparators
      bool operator == (HitState const& other) const { return _state == other._state; }
      bool operator != (HitState const& other) const { return !operator==(other); }
      bool operator < (HitState const& other ) const { return _state < other._state; }
      TSHState _state; //
    };
    // vector of HitStates
    typedef std::vector<HitState> HSV;
    // the state of all hits in a panel defines the panel state
    typedef std::vector<HitState> PanelState;
    // a vector to expliclty define all possible panel states
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
      XYZVectorF _udir;  // u direction in detector coordinates
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

    // utility function for combinatorics
    unsigned ipow(unsigned base, unsigned exp);
    bool setHitStates(PanelState const& pstate, TrkStrawHitVector& hits);
  }
}
#endif
