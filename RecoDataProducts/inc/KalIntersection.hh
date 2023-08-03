//
//  Class recording the intersection (geometric and physical) of a KinKal fit result with
//  a geometric surface.
//
#ifndef RecoDataProducts_KalIntersection_HH
#define RecoDataProducts_KalIntersection_HH
#include "Offline/DataProducts/inc/GenVector.hh"
#include "KinKal/General/ParticleStateEstimate.hh"
#include "Offline/KinKalGeom/inc/SurfaceId.hh"
#include <vector>
namespace mu2e {
  struct KalIntersection {
    XYZVectorF norm_; // surface unit normal at intersection point
    XYZVectorF bnom_; // Bfield at this intersection, needed to reconstitute trajectory
    KinKal::ParticleStateEstimate _pstate; // particle state at intersection point/time
    SurfaceId surface_; // which surface in the reco geometry was interestected
  };
  using KalIntersectionCollection = std::vector<KalIntersection>;
}

#endif
