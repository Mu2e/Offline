// predicate to extrapolate to the next intersection with the IPA
// track changes direction.
#ifndef Mu2eKinKal_ExtrapolateIPA_hh
#define Mu2eKinKal_ExtrapolateIPA_hh
#include "Offline/Mu2eKinKal/inc/KKTrack.hh"
#include "KinKal/General/TimeDir.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/Geometry/Intersection.hh"
#include "KinKal/Geometry/ParticleTrajectoryIntersect.hh"
#include <limits>
#include "cetlib_except/exception.h"
namespace mu2e {
  using KinKal::TimeDir;
  using KinKal::TimeRange;
  using KinKal::Intersection;
  class ExtrapolateIPA {
    public:
      using CylPtr = std::shared_ptr<KinKal::Cylinder>;
      ExtrapolateIPA() : maxDt_(-1.0), tol_(1e10) {}
      ExtrapolateIPA(double maxdt, double tol,CylPtr const& ipa) : maxDt_(maxdt), tol_(tol), ipa_(ipa) {}
      // interface for extrapolation
      double maxDt() const { return maxDt_; } // maximum time to extend the track
      double tolerance() const { return tol_; } // intersection tolerance
      CylPtr const& IPACylinder() const { return ipa_; }
      auto const& intersection() const { return inter_; }
      // extrapolation predicate: the track will be extrapolated until this predicate returns false, subject to the maximum time
      template <class KTRAJ> bool needsExtrapolation(KinKal::Track<KTRAJ> const& ktrk, TimeDir tdir, double time) const;
    private:
      double maxDt_; // maximum extrapolation time
      double tol_; // intersection tolerance (mm)
      CylPtr ipa_; // IPA cylinder
      mutable Intersection inter_; // most recent intersection
  };

  template <class KTRAJ> bool ExtrapolateIPA::needsExtrapolation(KinKal::Track<KTRAJ> const& ktrk, KinKal::TimeDir tdir, double time) const {
    // sanity check on initial position
    double dz = ktrk.fitTraj().position3(time).Z() - ipa_->center().Z();
    double halflen = ipa_->halfLength();
//    if(dz < -halflen || dz > halflen) return false; // need to qualify his in the case of reflection TODO
    // estimate time range based on local piece Z speed transit time
    double tz = 1.0/std::max(fabs(ktrk.fitTraj().velocity(time).Z())/(2*halflen),1.0/maxDt_); // protect against reflection (zero z speed)
    TimeRange trange = tdir == TimeDir::forwards ? TimeRange(time,time+tz) : TimeRange(time-tz,time);
    // update intersection
    std::cout << "IPA intersection range " << trange << std::endl;
    inter_ = KinKal::intersect(ktrk.fitTraj(),*ipa_,trange,tol_,tdir);
    std::cout << "IPA extrap start time " << time << " " << inter_.time_ << " " << inter_.onsurface_ << " " << inter_.inbounds_ << std::endl;
    return inter_.onsurface_ && inter_.inbounds_ && inter_.time_ < time;
  }
}
#endif
