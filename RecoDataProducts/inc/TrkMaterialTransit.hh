//
//  Estimate of the material traversed by a track fit (including its extrapolation) through a volume containing
//  distributed material, such as the stopping target.  The fit trajectory is integrated in time through the volume
//  bounding the material.  At each time, the probability that the true position is inside the volume is computed from the
//  trajectory position covariance projected onto the plane transverse to the velocity (along-track errors only re-label time).
//  The time integral of speed times that probability is the expected path length in the volume; multiplying by the material
//  fill factor (material volume / bounding volume) gives the expected path length through material.
//  Each contiguous sweep of the trajectory through the volume is recorded as a separate arc transit.
//
#ifndef RecoDataProducts_TrkMaterialTransit_hh
#define RecoDataProducts_TrkMaterialTransit_hh
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include <vector>
namespace mu2e {
  //
  //  TrkArcTransit: the transit estimate of one arc of the trajectory.  The payload is PROBABILISTIC: an arc is a contiguous time
  //  range over which the probability that the true trajectory is inside the bounding volume (given the fit covariance) exceeds a
  //  small threshold.  The nominal (central) trajectory may never actually cross into the bounding volume, nor into any of the
  //  material volumes (ie foils) inside it; such near-miss arcs have small but non-zero expected path lengths.  Conversely the
  //  material path length is an expectation averaged over the material fill factor, not the path through any particular piece of
  //  material.  The entry and exit times and surfaces describe where that probability rises above and falls below the threshold.
  //
  struct TrkArcTransit {
    double tentry_ = 0.0; // time the probability of being inside the volume rises above threshold (ns)
    double texit_ = 0.0; // time the probability of being inside the volume falls below threshold (ns)
    SurfaceId sentry_; // bounding surface limiting the probability at tentry_ (unknown if the trajectory starts inside the volume); the trajectory may not physically cross it
    SurfaceId sexit_; // bounding surface limiting the probability at texit_ (unknown if the trajectory ends inside the volume); the trajectory may not physically cross it
    float vollen_ = 0.0; // expected (probability-weighted) path length inside the bounding volume (mm)
    float matlen_ = 0.0; // expected path length through material = fill factor * vollen_ (mm)
    float mom_ = 0.0; // scalar momentum at the arc entrance (MeV/c)
    float eloss_ = 0.0; // energy loss estimated for the expected material path length matlen_ (MeV, positive for a loss)
    float elossrms_ = 0.0; // RMS of that energy loss (MeV)
    TrkArcTransit() {}
    double duration() const { return texit_ - tentry_; }
  };
  struct TrkMaterialTransit {
    SurfaceId sid_; // volume containing the material (ie ST_Foils)
    std::vector<TrkArcTransit> arcs_; // arc transits, in time order
    TrkMaterialTransit() {}
    TrkMaterialTransit(SurfaceId const& sid) : sid_(sid) {}
    auto const& arcs() const { return arcs_; }
    float matlen() const { float retval(0.0); for(auto const& arc : arcs_) retval += arc.matlen_; return retval; }
    float eloss() const { float retval(0.0); for(auto const& arc : arcs_) retval += arc.eloss_; return retval; }
  };
  using TrkMaterialTransitCollection = std::vector<TrkMaterialTransit>;
}
#endif
