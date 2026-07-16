// Predicate to extrapolate to the next intersection with passive planar material shells.
#ifndef Mu2eKinKal_ExtrapolatePlanes_hh
#define Mu2eKinKal_ExtrapolatePlanes_hh

#include "KinKal/General/TimeDir.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/Geometry/Intersection.hh"
#include "KinKal/Geometry/ParticleTrajectoryIntersect.hh"
#include "Offline/KinKalGeom/inc/KinKalGeom.hh"
#include "Offline/Mu2eKinKal/inc/ExtrapolateSurface.hh"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace mu2e {
  using KinKal::Intersection;
  using KinKal::TimeDir;
  using KinKal::TimeRange;
  using KinKal::timeDirSign;

  struct PlaneIntersection {
    using MaterialPlane = KinKalGeom::PassiveMaterialPlane;
    Intersection inter_;
    MaterialPlane const* plane_ = nullptr;
    PlaneIntersection() {}
    PlaneIntersection(Intersection const& inter, MaterialPlane const& plane) :
        inter_(inter), plane_(&plane) {}
  };

  class ExtrapolatePlanes : public ExtrapolateSurface<PlaneIntersection> {
    public:
      using Base = ExtrapolateSurface<PlaneIntersection>;
      using MaterialPlane = KinKalGeom::PassiveMaterialPlane;
      using MaterialPlaneCollection = KinKalGeom::PassiveMaterialPlaneCollection;

      ExtrapolatePlanes(double maxdt, double maxdtstep, double dptol, double intertol,
          double minv, MaterialPlaneCollection const& planes, int debug=0) :
        Base(maxdt,maxdtstep,dptol,intertol,minv,debug), planes_(planes) {}

      template <class KTRAJ> bool needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const;

    private:
      MaterialPlaneCollection const& planes_;
  };

  template <class KTRAJ> bool ExtrapolatePlanes::needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const {
    reset();
    if(planes_.empty()) return false;
    auto const& ktraj = tdir == TimeDir::forwards ? fittraj.back() : fittraj.front();
    static const double epsilon(1.0e-6);
    if(ktraj.range().range() <= epsilon) return true;
    auto stime = tdir == TimeDir::forwards ? ktraj.range().begin()+epsilon : ktraj.range().end()-epsilon;
    auto etime = tdir == TimeDir::forwards ? ktraj.range().end() : ktraj.range().begin();
    TimeRange trange(stime,etime,false);
    auto pos = ktraj.position3(stime);
    auto vel = ktraj.velocity(stime);
    double piece_span = std::fabs(stime - etime);
    bool keepGoing = false;
    for(auto const& plane : planes_) {
      double normvel = vel.Dot(plane.surface_->normal())*timeDirSign(tdir);
      double pdist = (plane.surface_->center()-pos).Dot(plane.surface_->normal());
      if(debug_ > 4)std::cout << "Passive plane extrap " << plane.sid_ << " normvel " << normvel << " pdist " << pdist << std::endl;
      if(std::fabs(normvel) < minvnorm_) continue;
      if(pdist*normvel < 0.0) continue;
      auto newinter = KinKal::intersect(fittraj,*plane.surface_,trange,intertol_,tdir);
      if(debug_ > 3)std::cout << plane.sid_ << " " << newinter << std::endl;
      if(newinter.good()) {
        inters_.emplace_back(newinter,plane);
      } else if(newinter.onsurface_ && newinter.inbounds_) {
        keepGoing |= trange.beyond(newinter.time_,tdir);
      } else if(!newinter.onsurface_ && tdir == TimeDir::backwards) {
        // a far plane beyond this (short) piece: the range-bounded stepIntersect can't reach it, so
        // use the analytic ray estimate (pdist/normvel) to keep marching outward (mirrors
        // ExtrapolateCRVRegion). This lets a tight KinematicLine-tail seed reach the concrete plane
        // instead of overshooting it (and the CRV) to maxDt.
        if(pdist/normvel > piece_span) keepGoing = true;
      }
    }
    SortIntersections isort(tdir);
    std::sort(inters_.begin(),inters_.end(),isort);
    // stop as soon as a crossing is cached so the KKExtrap loop records it HERE (no overshoot)
    if(!inters_.empty()) return false;
    return keepGoing;
  }
}

#endif
