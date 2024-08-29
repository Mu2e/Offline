// predicate for extrapolation to a given Z position in a given direction
// track changes direction.
#ifndef Mu2eKinKal_ExtrapolateToZ_hh
#define Mu2eKinKal_ExtrapolateToZ_hh
#include "Offline/Mu2eKinKal/inc/KKTrack.hh"
#include "KinKal/General/TimeDir.hh"
#include <limits>
#include "cetlib_except/exception.h"
namespace mu2e {
  using KinKal::TimeDir;
  class ExtrapolateToZ {
    public:
      ExtrapolateToZ() : maxDt_(-1.0), tol_(1e10), zval_(0.0) {}
      ExtrapolateToZ(double maxdt, double tol, double zval) : maxDt_(maxdt), tol_(tol), zval_(zval) {}
      // interface for extrapolation
      double maxDt() const { return maxDt_; } // maximum time to extend the track, WRT the time of the first(last) measurement
      double tolerance() const { return tol_; } // tolerance on fractional momentum change
      double zVal() const { return zval_; }
      // extrapolation predicate: the track will be extrapolated until this predicate returns false, subject to the maximum time
      template <class KTRAJ> bool needsExtrapolation(KinKal::Track<KTRAJ> const& ktrk, TimeDir tdir, double time) const;
    private:
      double maxDt_; // maximum extrapolation time
      double tol_; // momentum tolerance in BField domain
      double zval_; // z value targeted
  };

  template <class KTRAJ> bool ExtrapolateToZ::needsExtrapolation(KinKal::Track<KTRAJ> const& ktrk, KinKal::TimeDir tdir, double time) const {
    auto dir = ktrk.fitTraj().direction(time);
    auto dir0 = ktrk.fitTraj().direction(ktrk.fitTraj().t0());
    // check for reflection
    if(dir0.Z()*dir.Z() < 0.0) return false;
    auto pos = ktrk.fitTraj().position3(time);
    double zval = pos.Z();
    // sign the z velocity by the extrapolation direction
    double extrapdir = tdir == KinKal::TimeDir::forwards ? dir.Z() : -dir.Z();
    // stop when we get beyond the target value in Z
    if(extrapdir < 0){ // backwards extrapolation of downstream-going track, or forwards extrapolation of upstream-going track
      return zval > zval_;
    } else { // opposite
      return zval < zval_;
    }
  }
}
#endif
