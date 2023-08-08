//
//  Class recording the intersection (geometric and physical) of a KinKal fit result with
//  a geometric surface.
//
#ifndef RecoDataProducts_KalIntersection_HH
#define RecoDataProducts_KalIntersection_HH
#include "Offline/DataProducts/inc/GenVector.hh"
#include "KinKal/General/ParticleStateEstimate.hh"
#include "KinKal/Geometry/Intersection.hh"
#include "Offline/KinKalGeom/inc/SurfaceId.hh"
#include <vector>
namespace mu2e {
  struct KalIntersection {
    KinKal::ParticleStateEstimate pstate_; // particle state at intersection point/time
    XYZVectorF bnom_; // Bfield at this intersection, needed to reconstitute trajectory
    SurfaceId surfid_; // which surface in the reco geometry was interestected
    KinKal::IntersectFlag flag_; // intersection flag
    XYZVectorF norm_; // surface unit normal at intersection point
    KalIntersection(){}
    KalIntersection(KinKal::ParticleStateEstimate const& pstate, XYZVectorF const& bnom, SurfaceId const& surfid, KinKal::Intersection const& idata) : pstate_(pstate), bnom_(bnom), surfid_(surfid), flag_(idata.flag_), norm_(XYZVectorF(idata.norm_)) {}
// simple accessors
    double time() const { return pstate_.time(); }
    double mom() const { return pstate_.momentum(); }
    double momerr() const { return sqrt(pstate_.momentumVariance()); }
    KinKal::VEC3 momentum3() const { return pstate_.momentum3(); }
    KinKal::VEC3 velocity() const { return pstate_.velocity(); }
    KinKal::VEC3 position3() const { return pstate_.position3(); }
     // convert content to a LoopHelix
    KinKal::LoopHelix loopHelix() const { return KinKal::LoopHelix(pstate_, KinKal::VEC3(bnom_)); }
    // convert to a CentralHelix
    KinKal::CentralHelix centralHelix() const { return KinKal::CentralHelix(pstate_, KinKal::VEC3(bnom_)); }
    // convert to a KinematicLine
    KinKal::KinematicLine kinematicLine() const { return KinKal::KinematicLine(pstate_, KinKal::VEC3(bnom_)); }
 };
  using KalIntersectionCollection = std::vector<KalIntersection>;
}

#endif
