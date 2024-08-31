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
      ExtrapolateIPA() : maxDt_(-1.0), tol_(1e10),
      ipazmin_(std::numeric_limits<float>::max()),
      ipazmax_(-std::numeric_limits<float>::max()),debug_(0) {}

      ExtrapolateIPA(double maxdt, double tol,CylPtr const& ipa, int debug=0) : maxDt_(maxdt), tol_(tol), ipa_(ipa),
        ipazmin_( (ipa_->center() - ipa_->axis()*ipa_->halfLength()).Z()),
        ipazmax_( (ipa_->center() + ipa_->axis()*ipa_->halfLength()).Z()), debug_(debug) {}
      // interface for extrapolation
      double maxDt() const { return maxDt_; } // maximum time to extend the track
      double tolerance() const { return tol_; } // intersection tolerance
      CylPtr const& IPACylinder() const { return ipa_; }
      auto const& intersection() const { return inter_; }
      int debug() const { return debug_; }
      // extrapolation predicate: the track will be extrapolated until this predicate returns false, subject to the maximum time
      template <class KTRAJ> bool needsExtrapolation(KinKal::Track<KTRAJ> const& ktrk, TimeDir tdir, double time) const;
      // reset between tracks
      void reset() const { inter_ = Intersection(); }
    private:
      double maxDt_; // maximum extrapolation time
      double tol_; // intersection tolerance (mm)
      CylPtr ipa_; // IPA cylinder
      mutable Intersection inter_; // cache of most recent intersection
      // cache of IPA front and back Z positions
      double ipazmin_, ipazmax_;
      int debug_; // debug level
  };

  template <class KTRAJ> bool ExtrapolateIPA::needsExtrapolation(KinKal::Track<KTRAJ> const& ktrk, TimeDir tdir, double time) const {
    if(debug_ > 1)std::cout << "IPA extrap start time " << time << std::endl;
    auto vel = ktrk.fitTraj().velocity(time);
    auto pos = ktrk.fitTraj().position3(time);
    double zvel = vel.Z()*timeDirSign(tdir); // sign by extrapolation direction
    double zpos = pos.Z();
    // stop if we're heading away from the IPA
    if( (zvel > 0 && zpos > ipazmax_ ) || (zvel < 0 && zpos < ipazmin_)){
      reset(); // clear any cache
      if(debug_ > 1)std::cout << "Heading away from IPA: done" << std::endl;
      return false;
    }
    // if we're going in the right direction but haven't yet reached the IPA just keep going
    if( (zvel > 0 && zpos < ipazmin_) || (zvel < 0 && zpos > ipazmax_) ){
      reset();
      if(debug_ > 1)std::cout << "Heading towards IPA, z " << zpos<< std::endl;
      return true;
    }
    // if we get to here we need to test for an intersection with the actual cylinder
    // first, estimate the time range based on local piece Z speed transit time
    // Buffer by the range of the last piece to avoid missed edges.
    auto const& ktraj = tdir == TimeDir::forwards ? ktrk.fitTraj().back() : ktrk.fitTraj().front();
    static const double epsilon(1e-8);
    double dt = ktraj.range().range() - epsilon; // small difference to avoid re-intersecting
    double halflen = ipa_->halfLength();
    double tz = 1.0/std::max(fabs(zvel)/(2*halflen),1.0/maxDt_); // protect against reflection (zero z speed)
    TimeRange trange = tdir == TimeDir::forwards ? TimeRange(time-dt,time+tz) : TimeRange(time-tz,time+dt);
    // update intersection
    if(debug_ > 1)std::cout << "IPA intersection " << trange << std::endl;
    auto newinter = KinKal::intersect(ktrk.fitTraj(),*ipa_,trange,tol_,tdir);
    if(debug_ > 1)std::cout << "IPA extrap inter " << newinter.time_ << " " << newinter.onsurface_ << " " << newinter.inbounds_ << std::endl;
    bool goodextrap = newinter.onsurface_ && newinter.inbounds_;
    if(goodextrap){
      // if the cached intersection is valid, test this intersection time against it, and
      // if the new intersection time is the same as the last, keep extrapolating
      if(inter_.onsurface_ && inter_.inbounds_ && ( (tdir == TimeDir::forwards && newinter.time_ <= inter_.time_) ||
          (tdir == TimeDir::backwards && newinter.time_ >= inter_.time_) ) ) {
        if(debug_ > 1)std::cout << "Skipping duplicate intersection " << std::endl;
        return true;
      }
      // otherwise test if the trajectory extends to the intersection time yet or not. If so we are done
      if ( (tdir == TimeDir::forwards && newinter.time_ < time) ||
          (tdir == TimeDir::backwards && newinter.time_ > time ) ) {
        // update the cache
        inter_ = newinter;
        if(debug_ > 0)std::cout << "Good intersection found in range, time " << inter_.time_ << " z  " << inter_.pos_.Z() << std::endl;
        return false;
      }
      return true; // otherwise continue extrapolating
    } else {
      // no more intersections: keep extending in Z till we clear the IPA
      reset();
      if(debug_ > 1)std::cout << "Extrapolating to IPA edge, z " << zpos << std::endl;
      if(zvel > 0.0)
        return zpos < ipazmax_;
      else
        return zpos > ipazmin_;
    }
  }
}
#endif
