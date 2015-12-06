//
// class to define a set of allowed panel states and allow iterating through them
// This also interacts with the hits themselves
// Dave Brown
#ifndef mu2e_PanelAmbig_PanelStateIterator_HH
#define mu2e_PanelAmbig_PanelStateIterator_HH
#include "TrkReco/inc/PanelAmbigStructs.hh"

namespace mu2e {
  namespace PanelAmbig {

    class PanelStateIterator {
      public:
	PanelStateIterator() = default;
	// construct from a vector of TrkStrawHits and the list of allowed hit states.
	PanelStateIterator(TSHUIV const& uinfo, HSV const& allowed);
	// copy constructor
	PanelStateIterator(PanelStateIterator const& other) = default;
	PanelStateIterator& operator =(PanelStateIterator const& other) = default;
	// current state
	PanelState current() { return *_current; }
	// total # of states
	size_t nStates() const { return _allowedPS.size(); }
	// operate on the current state
	bool increment() { ++_current; return _current != _allowedPS.end(); }
	void reset() { _current = _allowedPS.begin(); }
      private:
	// helper functions;
	bool increment(HSV& hsv);
	bool increment(HitState& hs);
	void reset(HitState& hs);
	PSV::iterator _current; // current panel state
	TSHUIV _uinfo; // information about the panel hits
	PSV _allowedPS; // all allowed states for this panel
	HSV _allowedHS; // allowed hit states
    };
 
  } // PanelAmbig namespace
} // mu2e namespace
#endif
