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
    Int_t _ndem; // number of downstreameMinus tracks 
    Int_t _nuem; // number of upstreameMinus tracks 
    Int_t _ndmm; // number of downstreammuMinus tracks 
    Int_t _ndemc; // Number of calo clusters matched to the best dem track.
    Int_t _ndemo; // number of shared hits between primary and next-best track
    static std::string const& leafnames() { 
      static const std::string leaves =
	std::string("ndem/I:nuem/I:ndmm/I:ndemc/I:ndemo/I");
	return leaves;
    }
    void reset() {
      _ndem = _nuem = _ndmm = _ndemc = _ndemo = 0;
    }
  };
}
#endif
