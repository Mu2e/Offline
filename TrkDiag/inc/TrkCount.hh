#ifndef TrkCount_HH
#define TrkCount_HH
// root 
#include "Rtypes.h"
//
// Struct to count tracks and track-related quantities in an event
// Dave Brown, LBNL 7/8/2016
namespace mu2e 
{
  struct TrkCount {
    Int_t _nde; // number of downstreameMinus tracks 
    Int_t _nue; // number of upstreameMinus tracks 
    Int_t _ndm; // number of downstreammuMinus tracks 
    Int_t _ndec; // Number of calo clusters matched to the best dem track.
    Int_t _ndeo; // number of shared hits between primary and next-best track
    Int_t _ndmo; // number of shared hits between primary and muon-fit track
    static std::string const& leafnames() { 
      static const std::string leaves =
	std::string("nde/I:nue/I:ndmm/I:ndec/I:ndeo/I:ndmmo/I");
	return leaves;
    }
    void reset() {
      _nde = _nue = _ndm = _ndec = _ndeo = _ndmo = 0;
    }
  };
}
#endif
