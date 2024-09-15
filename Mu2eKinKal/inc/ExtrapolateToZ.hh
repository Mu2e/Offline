// predicate for extrapolation to a given Z position in a given direction
// track changes direction.
#ifndef Mu2eKinKal_ExtrapolateToZ_hh
#define Mu2eKinKal_ExtrapolateToZ_hh
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/General/TimeDir.hh"
#include <limits>
#include "cetlib_except/exception.h"
namespace mu2e {
  using KinKal::TimeDir;
  using KinKal::timeDirSign;
  class ExtrapolateToZ {
    public:
      ExtrapolateToZ() : maxDt_(-1.0), tol_(1e10), zval_(0.0), debug_(0) {}
      ExtrapolateToZ(double maxdt, double tol, double zval,int debug) : maxDt_(maxdt), tol_(tol), zval_(zval), debug_(debug) {}
      // interface for extrapolation
      double maxDt() const { return maxDt_; } // maximum time to extend the track, WRT the time of the first(last) measurement
      double tolerance() const { return tol_; } // tolerance on fractional momentum change
      double zVal() const { return zval_; }
      int debug() const { return debug_; }
      // extrapolation predicate: the track will be extrapolated until this predicate returns false, subject to the maximum time
      template <class KTRAJ> bool needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const;
    private:
      double maxDt_; // maximum extrapolation time
      double tol_; // momentum tolerance in BField domain
      double zval_; // z value targeted
      int debug_; // debug level
  };

  template <class KTRAJ> bool ExtrapolateToZ::needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, KinKal::TimeDir tdir) const {
    auto const& ktraj = tdir == TimeDir::forwards ? fittraj.back() : fittraj.front();
    auto time = tdir == TimeDir::forwards ? ktraj.range().end() : ktraj.range().begin();
    auto vel = ktraj.velocity(time);
    auto pos = ktraj.position3(time);
    double zvel = vel.Z()*timeDirSign(tdir); // sign by extrapolation direction
    double zval = pos.Z();
    auto const& bnom = ktraj.bnom(time);
    double zref = vel.R()*fabs(sin(bnom.Theta()));
    if(debug_ > 2)std::cout << "Z extrap start time " << time << " z " << zval << " zvel " << zvel << " zref " <<  zref << std::endl;
    // if z velocity is unreliable, continue
    if(fabs(zvel) < zref) return true;
    // stop if we're heading away from the target z
    if((zvel > 0 && zval > zval_ ) || (zvel < 0 && zval < zval_))return false;
    // stop when we get beyond the target value in Z
    if(zvel< 0){ // backwards extrapolation of downstream-going track, or forwards extrapolation of upstream-going track
      return zval > zval_;
    } else { // opposite
      return zval < zval_;
    }
  }
}
#endif
