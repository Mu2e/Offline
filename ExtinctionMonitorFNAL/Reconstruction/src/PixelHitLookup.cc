// Simulate pulse shape in the pixel time over threshold (ToT) circuit.
//
// Andrei Gaponenko, 2012

#include "ExtinctionMonitorFNAL/Reconstruction/inc/PixelHitLookup.hh"

#include <utility>

namespace mu2e {

  PixelHitLookup::PixelHitLookup(const ExtMonFNALRawHitCollection& hits)
    : origHits_(hits)
  {
    for(HitIndex i=0; i<hits.size(); ++i) {
      hitmap_.insert(std::make_pair(hits[i].pixelId(), i));
    }
  }

  ExtMonFNALRawHitCollection::size_type PixelHitLookup::findHit(const ExtMonFNALPixelId& id, int clock) const {
    typedef HitMap::const_iterator Iter;
    std::pair<Iter,Iter> range = hitmap_.equal_range(id);
    for(Iter i = range.first; i != range.second; ++i) {
      if(origHits_[i->second].clock() == clock) {
        // There can be only one hit for a given pixel at a given time,
        // don't have to look through the rest of the hits.
        return i->second;
      }
    }
    return HitIndex(-1);
  }

} // namespace mu2e
