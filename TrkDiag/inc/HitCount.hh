#ifndef HitCount_HH
#define HitCount_HH
// root 
#include "Rtypes.h"
//
// Struct to count hits in an event
// Dave Brown, LBNL 7/8/2016
namespace mu2e 
{
  struct HitCount {
    Int_t _nsd; // Total number of StrawDigis
    Int_t _nesel; // number of StrawHits flaged as passing energy deposition cuts
    Int_t _nrsel; // number of StrawHits flaged as passing radius cuts
    Int_t _ntsel; // number of StrawHits flaged as passing  time cuts
    Int_t _nbkg; // number of StrawHits flaged as produced by low-energy electrons (Compton electrons)
    Int_t _ntpk; // number of StrawHits assigned to a time peak
    Int_t _ncd, _ncrvd; // # calo digis, clusters, and crv digis
    static std::string const& leafnames() { 
      static const std::string leaves =
	std::string("nsh/I:nesel/I:nrsel/I:ntsel/I:nbkg/I:ntpk/I:ncd/I:ncrvd/I");
	return leaves;
    }
    void reset() {
      _nsd = _nesel = _nrsel = _ntsel = _nbkg = _ntpk = _ncd = _ncrvd = 0;
    }
  };
}
#endif
