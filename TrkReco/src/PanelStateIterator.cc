//
// class to define a set of allowed panel states and allow iterating through them
// This also interacts with the hits themselves
// Dave Brown

#include "TrkReco/inc/PanelStateIterator.hh"
#include <algorithm>

namespace mu2e {
  namespace PanelAmbig {
    
    PanelStateIterator::PanelStateIterator(TSHUIV const& uinfo, HSV const& allowed) : _uinfo(uinfo), _allowedHS(allowed) {
      // reserve enough space
      _allowedPS.reserve(ipow(_allowedHS.size(),_uinfo.size()));
      // setup an initial panel state.  The order follows the uinfo
      PanelState pstate;
      pstate.reserve(_uinfo.size());
      for(auto tshui : _uinfo) {
	if(tshui._use == TSHUInfo::free)
	// free hits are initialized to the first allowed state
	  pstate.push_back(_allowedHS.front());
	else
	// fixed or unused hits are kept at the initial state
	  pstate.push_back(tshui._hstate);
      }
      // this is the 1st panel state
      _allowedPS.push_back(pstate);
      // iterate over all possible states and record them
      while(increment(pstate)) {
	_allowedPS.push_back(pstate);
      }
      // set current to the 1st allowed state
      reset();
    }

    bool PanelStateIterator::increment(PanelState& pstate) {
      bool retval(false); // default is failure
      size_t ihit=0;
      do {
	if(_uinfo[ihit]._use == TSHUInfo::free){
	  HitState& hs = pstate[ihit];
	  if(increment(hs)){
	    retval = true;
	    break;
	  } else {
	  // we're at the end of allowed states for this hit; reset it and try incrementing the next
	    reset(hs);
	    ++ihit;
	  }
	} else
	  ++ihit; // skip hits not free to change
      } while(ihit < pstate.size());
      return retval;
    }

    bool PanelStateIterator::increment(HitState& hs) {
      bool retval(false);
      // find iterator to this state in the allowed states
      HSV::const_iterator ihs = std::find(_allowedHS.begin(), _allowedHS.end(),hs);
      // try incrementing: if successful, update the hit state, and we're done
      if(ihs != _allowedHS.end()) ++ihs;
      if(ihs != _allowedHS.end()){
	// success: update state
	hs = *ihs;
	retval = true;
      }
      return retval;
    }

    void PanelStateIterator::reset(HitState& hs) {
      hs = _allowedHS[0];
    }

  } // PanelAmbig namespace
} // mu2e namespace

