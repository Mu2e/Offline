//
// simple struct to record event-level information
//
#ifndef EventInfo_HH
#define EventInfo_HH
#include "Rtypes.h"
#include <string>
namespace mu2e
{
  struct EventInfo {
    Int_t _eventid, _runid, _subrunid; // run/event identification
    Float_t _evtwt, _beamwt, _genwt; // event weights: total, beam, generator
    Int_t _nprotons; // # of protons on target for this microbunch
    static std::string const& leafnames() { 
      static const std::string leaves =
	std::string("eventid/I:runid/I:subrunid/I:") +
	std::string("evtwt/F:beamwt/F:genwt/F:") +
	std::string("nprotons/I");
	return leaves;
    }
    void reset() {
      _eventid = _runid = _subrunid = _nprotons = 0;
      _evtwt = _beamwt = _genwt = 1.0;
    }
  };
}
#endif
