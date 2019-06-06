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
    static const int MAX_COUNTS = 50;
    //    Int_t _ncounts;
    Int_t _counts[MAX_COUNTS]; // number of tracks in collection
    //    Int_t _overlaps[MAX_COUNTS]; // number of shared hits between candidate track and this track

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
    const std::string leafnames(const std::vector<std::string>& trkbranches) { 
      std::string leaves;// = "ncnt/I:";
      for (std::vector<std::string>::const_iterator i_trkbranch = trkbranches.begin(); i_trkbranch != trkbranches.end(); ++i_trkbranch) {
	leaves += "n" + *i_trkbranch + "/I";
	//      }
	//      for (std::vector<std::string>::const_iterator i_trkbranch = trkbranches.begin(); i_trkbranch != trkbranches.end(); ++i_trkbranch) {
	//	leaves += "n" + *i_trkbranch + "o/I";
	if (i_trkbranch != trkbranches.end()-1) {
	  leaves += ":";
	}
      }
      return leaves;
    }
    void reset() {
      for (auto& i_count : _counts) {
	i_count = 0;
      }
      //      for (auto& i_overlap : _overlaps) {
      //	i_overlap = 0;
      //      }
      _nde = _nue = _ndm = _ndec = _ndeo = _ndmo = 0;
    }
  };
}
#endif
