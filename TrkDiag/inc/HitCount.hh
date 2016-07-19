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
    Int_t _nsh; // Total number of StrawHits
    Int_t _nesel; // number of StrawHits flaged as passing energy deposition cuts
    Int_t _nrsel; // number of StrawHits flaged as passing radius cuts
    Int_t _ntsel; // number of StrawHits flaged as passing  time cuts
    Int_t _nbkg; // number of StrawHits flaged as produced by low-energy electrons (Compton electrons)
    Int_t _nster; // number of StrawHits with stereo position information
    Int_t _ntdiv; // number of StrawHits with time division position information
    Int_t _ntpk; // number of StrawHits assigned to a time peak
    Int_t _nxt; // number of StrawHits flaged as cross-talk
    static std::string const& leafnames() { 
      static const std::string leaves =
	std::string("nsh/I:nesel/I:nrsel/I:ntsel/I:nbkg/I:nster/I:ntdiv/I:ntpk/I:nxt/I");
	return leaves;
    }
    void reset() {
      _nsh = _nesel = _nrsel = _ntsel = _nxt = _nbkg = _nster = _ntdiv = _ntpk = 0;
    }
  };
}
#endif
