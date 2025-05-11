// predicate to extrapolate to the next intersection with the IPA
// track changes direction.
#ifndef Mu2eKinKal_ExtrapolateIPA_hh
#define Mu2eKinKal_ExtrapolateIPA_hh
#include "KinKal/Trajectory/ParticleTrajectory.hh"
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
      ExtrapolateIPA() : maxDt_(-1.0), dptol_(1e10), intertol_(1e10),
      zmin_(std::numeric_limits<float>::max()),
      zmax_(-std::numeric_limits<float>::max()),debug_(0) {}

      ExtrapolateIPA(double maxdt, double dptol,double intertol, CylPtr const& ipa, int debug=0) :
        maxDt_(maxdt), dptol_(dptol), intertol_(intertol), ipa_(ipa),
        zmin_( (ipa_->center() - ipa_->axis()*ipa_->halfLength()).Z()),
        zmax_( (ipa_->center() + ipa_->axis()*ipa_->halfLength()).Z()), debug_(debug) {}
      // interface for extrapolation
      double maxDt() const { return maxDt_; }
      double dpTolerance() const { return dptol_; }
      double interTolerance() const { return intertol_; }
      CylPtr const& IPACylinder() const { return ipa_; }
      auto const& intersection() const { return inter_; }
      double zmin() const { return zmin_; }
      double zmax() const { return zmax_; }
      int debug() const { return debug_; }
      // extrapolation predicate: the track will be extrapolated until this predicate returns false, subject to the maximum time
      template <class KTRAJ> bool needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const;
      // reset between tracks
      void reset() const { inter_ = Intersection(); }
    private:
      double maxDt_; // maximum extrapolation time
      double dptol_; // fractional momentum tolerance
      double intertol_; // intersection tolerance (mm)
      CylPtr ipa_; // IPA cylinder
      mutable Intersection inter_; // cache of most recent intersection
      // cache of IPA front and back Z positions
      double zmin_, zmax_;
      int debug_; // debug level
  };

  template <class KTRAJ> bool ExtrapolateIPA::needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const {
    // we are answering the question: did the segment last added to this extrapolated track hit the IPA or not?
    // if so, stop extrapolating (for now). If not, and if we're still inside or heading towards the IPA, keep going.
    auto const& ktraj = tdir == TimeDir::forwards ? fittraj.back() : fittraj.front();
    // add a small buffer to the test range to prevent re-intersection with the same piece
    static const double epsilon(1e-7); // small step to avoid re-intersecting
    if(ktraj.range().range() <= epsilon) return true; // keep going if the step is very small
    auto stime = tdir == TimeDir::forwards ? ktraj.range().begin()+epsilon : ktraj.range().end()-epsilon;
    auto etime = tdir == TimeDir::forwards ? ktraj.range().end() : ktraj.range().begin();
    auto vel = ktraj.speed(stime)*ktraj.axis(stime).direction();// use helix axis to define the velocity
    auto spos = ktraj.position3(stime);
    auto epos = ktraj.position3(etime);
    double zvel = vel.Z()*timeDirSign(tdir); // sign by extrapolation direction
    if(debug_ > 2)std::cout << "IPA extrap start time " << stime << " start z " << spos.Z() << " end z " << epos.Z() << " zvel " << zvel << std::endl;
    // stop if the particle is heading away from the IPA
    if( (zvel > 0 && spos.Z() > zmax_ ) || (zvel < 0 && spos.Z() < zmin_)){
      reset(); // clear any cache
      if(debug_ > 1)std::cout << "Heading away from IPA: done" << std::endl;
      return false;
    }
    // if the particle is going in the right direction but hasn't yet reached the IPA just keep going
    if( (zvel > 0 && epos.Z() < zmin_) || (zvel < 0 && epos.Z() > zmax_) ){
      reset();
      if(debug_ > 2)std::cout << "Heading towards IPA, z " << spos.Z()<< std::endl;
      return true;
    }
    // if we get to here we need to test for an intersection with the actual cylinder. Make sure the range is positive definite
    auto trange = tdir == TimeDir::forwards ? TimeRange(stime,etime) : TimeRange(etime,stime);
    Intersection newinter = KinKal::intersect(fittraj,*ipa_,trange,intertol_,tdir);
    if(debug_ > 2)std::cout << "IPA " << newinter << std::endl;
    if(newinter.good()){
      // update the cache
      inter_ = newinter;
      if(debug_ > 0)std::cout << "Good IPA " << newinter << std::endl;
      return false;
    } else {
      // no more intersections: keep extending in Z till we clear the IPA
      reset();
      if(debug_ > 1)std::cout << "Extrapolating to IPA edge, z " << spos.Z() << std::endl;
      if(zvel > 0.0)
        return spos.Z() < zmax_;
      else
        return spos.Z() > zmin_;
    }
  }
}
#endif
