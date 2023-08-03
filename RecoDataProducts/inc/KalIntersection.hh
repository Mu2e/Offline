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
    KinKal::ParticleStateEstimate pstate_; // particle state at intersection point/time
    XYZVectorF bnom_; // Bfield at this intersection, needed to reconstitute trajectory
    SurfaceId surfid_; // which surface in the reco geometry was interestected
    XYZVectorF norm_; // surface unit normal at intersection point
    KalIntersection(){}
    KalIntersection(KinKal::ParticleStateEstimate const& pstate, XYZVectorF const& bnom, SurfaceId const& surfid, XYZVectorF const& norm) : pstate_(pstate), bnom_(bnom), surfid_(surfid), norm_(norm) {}
  };
  using KalIntersectionCollection = std::vector<KalIntersection>;
}

#endif
