// predicate for extrapolation to a given Z position in a given direction. Stops if the
// track changes direction.
#ifndef Mu2eKinKal_ExtrapolateToZ_hh
#define Mu2eKinKal_ExtrapolateToZ_hh
#include "Offline/Mu2eKinKal/inc/KKTrack.hh"
#include "KinKal/General/TimeDir.hh"
#include <limits>
#include "cetlib_except/exception.h"
namespace mu2e {
  class ExtrapolateToZ {
    public:
      ExtrapolateToZ() : maxDt_(-1.0), tol_(1e10),
      zmin_(std::numeric_limits<double>::max()),
      zmax_(-std::numeric_limits<double>::max()){}
      ExtrapolateToZ(double maxdt, double tol, double zmin, double zmax) :
        maxDt_(maxdt), tol_(tol), zmin_(zmin), zmax_(zmax) {
          if(zmin >= zmax) throw cet::exception("RECO")<<"Mu2eKinKal::ExtrapolateToZ: range configuration error, zmin "<< zmin << " zmax " << zmax <<endl;
        }
      // interface for extrapolation
      double maxDt() const { return maxDt_; } // maximum time to extend the track, WRT the time of the first(last) measurement
      double tolerance() const { return tol_; } // tolerance on fractional momentum change
      // extrapolation predicate: the track will be extrapolated until this predicate returns false, subject to the maximum time
      template <class KTRAJ> bool needsExtrapolation(KinKal::Track<KTRAJ> const& ktrk, KinKal::TimeDir tdir, double time) const;
    private:
      double maxDt_; // maximum extrapolation time
      double tol_; // momentum tolerance in BField domain
      double zmin_; // minimum z value required
      double zmax_; // maximum z value required
  };

  template <class KTRAJ> bool ExtrapolateToZ::needsExtrapolation(KinKal::Track<KTRAJ> const& kktrk, KinKal::TimeDir tdir, double time) const {
    auto pos = kktrk.fitTraj().position3(time);
    auto dir = kktrk.fitTraj().direction(time);
    double zval = pos.Z();
    // sign the z velocity by the extrapolation direction
    double extrapdir = tdir == KinKal::TimeDir::forwards ? dir.Z() : -dir.Z();
    // stop when we get outside the range. This works with reflections
    if(extrapdir < 0){ // backwards extrapolation of downstream-going track, or forwards extrapolation of upstream-going track
      return zval > zmin_;
    } else {
      return zval > zmax_;
    }
  }
}
#endif
