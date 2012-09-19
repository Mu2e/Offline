// A class to look up whether a given pixel was hit at a given time.
//
// Andrei Gaponenko, 2012

#ifndef ExtinctionMonitorFNAL_Reconstruction_PixelHitLookup_hh
#define ExtinctionMonitorFNAL_Reconstruction_PixelHitLookup_hh

#include <map>

// can't forward declare a typedef
#include "RecoDataProducts/inc/ExtMonFNALRawHitCollection.hh"

namespace mu2e {
  class ExtMonFNALPixelId;
  class ExtMonFNALRawHit;

  class PixelHitLookup {
  public:

    explicit PixelHitLookup(const ExtMonFNALRawHitCollection& hits);

    // returns index of the hit in the original collection or (-1) if not hit.
    typedef ExtMonFNALRawHitCollection::size_type HitIndex;
    HitIndex findHit(const ExtMonFNALPixelId& id, int clock) const;

  private:
    typedef std::multimap<ExtMonFNALPixelId, HitIndex> HitMap;
    HitMap hitmap_;
    const ExtMonFNALRawHitCollection &origHits_;
  };

} // namespace mu2e

#endif/*ExtinctionMonitorFNAL_Reconstruction_PixelHitLookup_hh*/
