// predicate for extrapolation to a given Z position in a given direction
// track changes direction.
#ifndef Mu2eKinKal_ExtrapolateToZ_hh
#define Mu2eKinKal_ExtrapolateToZ_hh
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/General/TimeDir.hh"
#include <limits>
#include <cmath>
#include "cetlib_except/exception.h"
namespace mu2e {
  using KinKal::TimeDir;
  using KinKal::timeDirSign;
  class ExtrapolateToZ {
    public:
      ExtrapolateToZ(double maxdt, double maxdtstep, double tol, double zval,int debug,
          double maxradius = std::numeric_limits<double>::max()) :
        maxDt_(maxdt), maxDtStep_(maxdtstep), dptol_(tol), zval_(zval), debug_(debug), maxradius_(maxradius) {}
      // interface for extrapolation
      double maxDt() const { return maxDt_; } // maximum time to extend the track, WRT the time of the first(last) measurement
      double maxDtStep() const { return maxDtStep_; }
      double dpTolerance() const { return dptol_; } // tolerance on fractional momentum change
      double zVal() const { return zval_; }
      int debug() const { return debug_; }
      // extrapolation predicate: the track will be extrapolated until this predicate returns false, subject to the maximum time
      template <class KTRAJ> bool needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const;
    private:
      double maxDt_ = -1; // maximum extrapolation time
      double maxDtStep_ = -1; // maximum extrapolation time step in a single iteration
      double dptol_ = 1e10; // fractional momentum tolerance in BField domain
      double zval_ = 0; // z value targeted
      int debug_ = 0; // debug level
      double maxradius_ = std::numeric_limits<double>::max(); // max transverse radius rho=sqrt(x^2+y^2) (mm)
                            // the track may reach; stops near-vertical (small-v_z) cosmics from chasing the
                            // constant-z target plane out to large rho. Default = no limit.
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
    double rho = std::sqrt(pos.X()*pos.X() + pos.Y()*pos.Y());
    if(debug_ > 2)std::cout << "Z extrap start time " << time << " z " << zval << " zvel " << zvel << " zref " <<  zref << " rho " << rho << std::endl;
    // stop if we have reached the maximum allowed transverse radius: a near-vertical (small-v_z) track
    // would otherwise keep chasing the constant-z target plane out to large rho (e.g. a cosmic running
    // metres in y), pre-extending the trajectory past the DS before the DS-material pass. Default = no limit.
    if(rho > maxradius_) return false;
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
