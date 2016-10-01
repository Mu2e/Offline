#ifndef RecoDataProducts_HelixHit_hh
#define RecoDataProducts_HelixHit_hh
//
// subclass of StrawHitPosition used for helix finding.  It just adds
// the resolved helix azimuth (relative to the helix axis) and a pointer
// to the original straw hit (since the collection will be sparse).
//
// Original author David Brown
//
// Mu2e includes
#include "RecoDataProducts/inc/StrawHitPosition.hh"
#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include <vector>

namespace mu2e {

  struct HelixHit : public StrawHitPosition {
    HelixHit() : _shidx(0), _phi(0.0) {}
    HelixHit(StrawHitPosition const& shp, StrawHitIndex shidx=0, Float_t phi=0.0) : StrawHitPosition(shp), _shidx(shidx), _phi(phi) {}
    StrawHitIndex _shidx; // index to the straw hit
    Float_t _phi; // resolved azimuth of this hit WRT the helix axis (circle center)
  };
  // define the collection type
  typedef std::vector<mu2e::HelixHit> HelixHitCollection;
}
#endif


