// predicate to extrapolate to the next intersection with the IPA
// track changes direction.
#ifndef Mu2eKinKal_ExtrapolateIPA_hh
#define Mu2eKinKal_ExtrapolateIPA_hh
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/General/TimeDir.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/Geometry/Intersection.hh"
#include "KinKal/Geometry/Cylinder.hh"
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

      ExtrapolateIPA(double maxdt, double maxdtstep, double dptol,double intertol, CylPtr const& ipa, int debug=0) :
        maxDt_(maxdt), maxDtStep_(maxdtstep), dptol_(dptol), intertol_(intertol), ipa_(ipa),
        zmin_( (ipa_->center() - ipa_->axis()*ipa_->halfLength()).Z()),
        zmax_( (ipa_->center() + ipa_->axis()*ipa_->halfLength()).Z()), debug_(debug) {}
      // interface for extrapolation
      double maxDt() const { return maxDt_; }
      double maxDtStep() const { return maxDtStep_; }
      double dpTolerance() const { return dptol_; }
      double interTolerance() const { return intertol_; }
      CylPtr const& IPACylinder() const { return ipa_; }
      auto const& intersection() const { return inter_; }
      double zmin() const { return zmin_; }
      double zmax() const { return zmax_; }
      double zmid() const { return 0.5*(zmin_+zmax_); }
      int debug() const { return debug_; }
      // extrapolation predicate: the track will be extrapolated until this predicate returns false, subject to the maximum time
      template <class KTRAJ> bool needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const;
      // reset between tracks
      void reset() const { inter_ = Intersection(); }
    private:
      double maxDt_ = -1; // maximum extrapolation time
      double maxDtStep_ = -1; // maximum extrapolation time step in a single iteration
      double dptol_ = 1e10; // fractional momentum tolerance
      double intertol_ = 1e10; // intersection tolerance (mm)
      CylPtr ipa_; // IPA cylinder
      mutable Intersection inter_; // cache of most recent intersection
      // cache of IPA front and back Z positions
      double zmin_ = std::numeric_limits<double>::max();
      double zmax_ = std::numeric_limits<double>::lowest();
      int debug_ = 0; // debug level
  };

  template <class KTRAJ> bool ExtrapolateIPA::needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const {
    // we are answering the question: did the segment last added to this extrapolated trajectory hit the IPA or not?
    // if so, stop extrapolating (for now). If not, and if we're still inside or heading towards the IPA, keep going, otherwise stop
    reset(); // clear any cache
    bool retval(true); // by default keep going
    auto const& ktraj = tdir == TimeDir::forwards ? fittraj.back() : fittraj.front(); // most recently added segment
    // add a small buffer to the test range to prevent re-intersection with the same piece
    static const double epsilon(1e-7); // small step to avoid re-intersecting
    if(ktraj.range().range() <= epsilon) return true; // keep going if the step is very small
    auto trange = tdir == TimeDir::forwards ?
      TimeRange(ktraj.range().begin()+epsilon, ktraj.range().end()) :
      TimeRange(ktraj.range().begin(), ktraj.range().end()-epsilon);
    auto spos = ktraj.position3(trange.begin());
    auto epos = ktraj.position3(trange.end());
    // if either end of the segment is inside the z range of the IPA, test for an in-range interestection
    if( (spos.Z() > zmin_ && spos.Z() < zmax_) || (epos.Z() > zmin_ && epos.Z() < zmax_) ){
      Intersection newinter = KinKal::intersect(ktraj,*ipa_,trange,intertol_,tdir);
      if(debug_ > 2)std::cout << "IPA " << newinter << std::endl;
      if(newinter.good()){
        // update the cache
        inter_ = newinter;
        if(debug_ > 0)std::cout << "Good IPA " << newinter << ", Stopping " << std::endl;
        retval = false;
      }
    } else {
      // otherwise, if the trajectory is heading towards the IPA keep going
      auto vel = ktraj.velocity(trange.end())*timeDirSign(tdir); // sign velocity by extrapolation direction
      double vb = vel.Dot(ktraj.bnom().Unit()); // project along the BField axis. This is insensitive to reflection.
      if( (epos.Z() > zmax_ && vb < 0) || (epos.Z() < zmin_ && vb > 0) ) {
        if(debug_ > 1)std::cout << "Extrapolating towards IPA, z " << epos.Z() << " bvel " << vb << std::endl;
      } else {
        retval = false;
        if(debug_ > 1)std::cout << "Heading away from IPA, z " << epos.Z() << " bvel " << vb << std::endl;
      }
    }
    return retval;
  }
}
#endif
