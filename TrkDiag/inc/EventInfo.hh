//
// simple struct to record event-level information
// Dave Brown (LBNL)
//
#ifndef EventInfo_HH
#define EventInfo_HH
#include "Rtypes.h"
#include <string>
namespace mu2e
{
  struct EventInfo {
    Int_t _eventid, _runid, _subrunid; // run/event identification
    Int_t _nprotons; // # of protons on target for this microbunch
    static std::string const& leafnames() { 
      static const std::string leaves =
	std::string("eventid/I:runid/I:subrunid/I:") +
	std::string("nprotons/I");
	return leaves;
    }
    void reset() {
      _eventid = _runid = _subrunid = _nprotons = 0;
    }
  };
}
#endif
