//
//  Class recording the intersection (geometric and physical) of a KinKal fit result with
//  a geometric surface.
//
#ifndef RecoDataProducts_KalIntersection_HH
#define RecoDataProducts_KalIntersection_HH
#include "Offline/DataProducts/inc/GenVector.hh"
#include "KinKal/General/ParticleStateEstimate.hh"
#include "KinKal/Geometry/Intersection.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include <vector>
namespace mu2e {
  struct KalIntersection {
    KinKal::ParticleStateEstimate pstate_; // particle state at intersection point/time
    XYZVectorF bnom_; // Bfield at this intersection, needed to reconstitute trajectory
    SurfaceId surfid_; // which surface in the reco geometry was interestected
    KinKal::Intersection kkinter_; // kinkal intersection
    double dP_ = 0.0; // estimated scalar momentum change in this intersection
    KalIntersection(){}
    KalIntersection(KinKal::ParticleStateEstimate const& pstate, XYZVectorF const& bnom, SurfaceId const& surfid, KinKal::Intersection const& kkinter,double dP=0.0) : pstate_(pstate), bnom_(bnom), surfid_(surfid), kkinter_(kkinter) ,dP_(dP){}
// simple accessors
    auto const& surfaceId() const { return surfid_; }
    double time() const { return pstate_.time(); }
    double mom() const { return pstate_.momentum(); }
    double momerr() const { return sqrt(pstate_.momentumVariance()); }
    double dMom() const { return dP_; }
    KinKal::VEC3 momentum3() const { return pstate_.momentum3(); }
    KinKal::VEC3 velocity() const { return pstate_.velocity(); }
    KinKal::VEC3 position3() const { return pstate_.position3(); }
    auto const & intersection() const { return kkinter_; }
    auto inBounds() const { return kkinter_.inbounds_; }
    auto const& surfaceNormal() const { return kkinter_.norm_; }
    auto gap() const { return kkinter_.gap_; }
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
