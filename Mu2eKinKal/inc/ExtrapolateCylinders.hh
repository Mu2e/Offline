// Predicate to extrapolate to the next intersection with passive cylindrical material shells.
#ifndef Mu2eKinKal_ExtrapolateCylinders_hh
#define Mu2eKinKal_ExtrapolateCylinders_hh

#include "KinKal/General/TimeDir.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/Geometry/Intersection.hh"
#include "KinKal/Geometry/ParticleTrajectoryIntersect.hh"
#include "Offline/KinKalGeom/inc/DetectorSolenoid.hh"
#include "Offline/Mu2eKinKal/inc/ExtrapolateSurface.hh"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace mu2e {
  using KinKal::Intersection;
  using KinKal::TimeDir;
  using KinKal::TimeRange;
  using KinKal::timeDirSign;

  struct CylinderIntersection {
    using MaterialCylinder = KKGeom::DetectorSolenoid::MaterialCylinder;
    Intersection inter_;
    MaterialCylinder const* cylinder_ = nullptr;
    CylinderIntersection() {}
    CylinderIntersection(Intersection const& inter, MaterialCylinder const& cylinder) :
        inter_(inter), cylinder_(&cylinder) {}
  };

  class ExtrapolateCylinders : public ExtrapolateSurface<CylinderIntersection> {
    public:
      using Base = ExtrapolateSurface<CylinderIntersection>;
      using MaterialCylinder = KKGeom::DetectorSolenoid::MaterialCylinder;
      using MaterialCylinderCollection = KKGeom::DetectorSolenoid::MaterialCylinderCollection;

      ExtrapolateCylinders(double maxdt, double maxdtstep, double dptol, double intertol,
          double minv, MaterialCylinderCollection const& cylinders, int debug=0) :
        Base(maxdt,maxdtstep,dptol,intertol,minv,debug), cylinders_(cylinders) {
        for(auto const& cylinder : cylinders_) {
          maxRadius_ = std::max(maxRadius_, cylinder.surface_->radius() + 0.5*cylinder.thickness_);
        }
      }

      template <class KTRAJ> bool needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const;

    private:
      double minShellCosine_ = 0.1; // min |cos(incidence)| to accept a shell crossing (~within 84 deg of normal);
                                    // below this the thin-shell path length blows up (near-tangent), so skip it
      MaterialCylinderCollection const& cylinders_;
      double maxRadius_ = 0.0;
  };

  template <class KTRAJ> bool ExtrapolateCylinders::needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const {
    reset();
    if(cylinders_.empty()) return false;
    auto const& ktraj = tdir == TimeDir::forwards ? fittraj.back() : fittraj.front();
    static const double epsilon(1.0e-6);
    if(ktraj.range().range() <= epsilon) return true;
    double stime = tdir == TimeDir::forwards ? ktraj.range().begin()+epsilon : ktraj.range().end()-epsilon;
    double etime = tdir == TimeDir::forwards ? ktraj.range().end() : ktraj.range().begin();
    // March the shell search across the end piece in BOUNDED sub-windows (capped by maxDtStep) from the
    // near end outward, rather than one intersect over the whole piece. The thin concentric shells are
    // closely spaced (~960-1315mm), so a single intersect over a piece that spans several of them can
    // step over the inner ones (they were under-found, innermost worst). Marching brackets the nearest
    // crossings and records them HERE, with no overstep; the KKExtrap loop banks/de-dups and re-extends.
    double const sgn = (etime >= stime) ? 1.0 : -1.0;
    double const dtwin = (maxDtStep_ > 0.0) ? maxDtStep_ : std::fabs(etime - stime);
    double t0 = stime;
    while(sgn*(etime - t0) > epsilon) {
      // stop once we are outside the shells AND moving further out: any shell is behind us
      auto p0 = ktraj.position3(t0);
      double const rho0 = std::sqrt(p0.X()*p0.X() + p0.Y()*p0.Y());
      auto v0 = ktraj.velocity(t0);
      double const vrho0 = rho0 > 0.0 ? (p0.X()*v0.X() + p0.Y()*v0.Y())/rho0 : 0.0;
      if(rho0 > maxRadius_ && vrho0*timeDirSign(tdir) > 0.0) return false;
      double const t1 = t0 + sgn*std::min(dtwin, std::fabs(etime - t0));
      TimeRange win = (sgn > 0.0) ? TimeRange(t0,t1,false) : TimeRange(t1,t0,false);
      for(auto const& cylinder : cylinders_) {
        auto newinter = KinKal::intersect(fittraj,*cylinder.surface_,win,intertol_,tdir);
        if(debug_ > 3)std::cout << cylinder.sid_ << " " << newinter << std::endl;
        if(newinter.good()) {
          double const normvel = std::fabs(newinter.norm_.Dot(newinter.pdir_));
          // thin-shell approximation: skip near-tangent crossings (unbounded path through a finite shell)
          if(normvel > std::max(minvnorm_,minShellCosine_)) inters_.emplace_back(newinter,cylinder);
          else if(debug_ > 1) std::cout << "Skipping grazing DS material intersection " << cylinder.sid_
            << " normvel " << normvel << std::endl;
        }
      }
      if(!inters_.empty()) {
        SortIntersections isort(tdir);
        std::sort(inters_.begin(),inters_.end(),isort);
        return false; // nearest shell crossing(s) found and recorded here -- no overstep
      }
      t0 = t1;
    }
    // reached the piece end with no crossing while still approaching the shells: extend and retry
    return true;
  }
}

#endif
