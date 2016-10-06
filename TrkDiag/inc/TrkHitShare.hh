#ifndef TrkHitShare_HH
#define TrkHitShare_HH
#include "Rtypes.h"
namespace mu2e
{
  // struct to look for hit sharing between tracks
  struct TrkHitShare {
    UInt_t _trk1; // index to primary track (== track with the most active hits)
    UInt_t _trk2; // index to secondary track (has shared hits)
    UInt_t _nhshared; // # of active hits shared between 2 tracks
    Float_t _f1; // fraction of shared active hits in primary track
    Float_t _f2; // fraction of shared active hits in secondary track
    static std::string leafnames() { static std::string leaves;
      leaves = std::string("trk1/i:trk2/i:nhshared/i:frac1/F:frac2/F");
      return leaves;
    }
  };
}
#endif

