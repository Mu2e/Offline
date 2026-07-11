// predicate for extrapolation to a given transverse radius rho = sqrt(x^2+y^2), bounded in z.
// This is the radial analog of ExtrapolateToZ: it targets a constant-rho (cylindrical) surface rather
// than a constant-z plane. That is the natural tracker exit for a CentralHelix (cosmic) track, which
// traverses the tracker radially (mostly in y) and leaves through the cylindrical side, NOT through the
// z-end planes. The [minz,maxz] bounds keep a near-horizontal (large v_z) cosmic from chasing the target
// radius far down the beamline, past the tracker ends -- the symmetric guard to ExtrapolateToZ's behaviour
// for near-vertical cosmics.
#ifndef Mu2eKinKal_ExtrapolateToR_hh
#define Mu2eKinKal_ExtrapolateToR_hh
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/General/TimeDir.hh"
#include <limits>
#include <cmath>
#include <iostream>
#include "cetlib_except/exception.h"
namespace mu2e {
  using KinKal::TimeDir;
  using KinKal::timeDirSign;
  class ExtrapolateToR {
    public:
      ExtrapolateToR(double maxdt, double maxdtstep, double tol, double rval, double minz, double maxz, int debug) :
        maxDt_(maxdt), maxDtStep_(maxdtstep), dptol_(tol), rval_(rval), minz_(minz), maxz_(maxz), debug_(debug) {}
      // interface for extrapolation
      double maxDt() const { return maxDt_; } // maximum time to extend the track, WRT the time of the first(last) measurement
      double maxDtStep() const { return maxDtStep_; }
      double dpTolerance() const { return dptol_; } // tolerance on fractional momentum change
      double rVal() const { return rval_; }
      double minZ() const { return minz_; }
      double maxZ() const { return maxz_; }
      int debug() const { return debug_; }
      // extrapolation predicate: the track will be extrapolated until this predicate returns false, subject to the maximum time
      template <class KTRAJ> bool needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const;
    private:
      double maxDt_ = -1; // maximum extrapolation time
      double maxDtStep_ = -1; // maximum extrapolation time step in a single iteration
      double dptol_ = 1e10; // fractional momentum tolerance in BField domain
      double rval_ = 0; // target transverse radius (mm)
      double minz_ = -std::numeric_limits<double>::max(); // minimum z allowed during extrapolation (mm)
      double maxz_ =  std::numeric_limits<double>::max(); // maximum z allowed during extrapolation (mm)
      int debug_ = 0; // debug level
  };

  template <class KTRAJ> bool ExtrapolateToR::needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, KinKal::TimeDir tdir) const {
    auto const& ktraj = tdir == TimeDir::forwards ? fittraj.back() : fittraj.front();
    auto time = tdir == TimeDir::forwards ? ktraj.range().end() : ktraj.range().begin();
    auto vel = ktraj.velocity(time);
    auto pos = ktraj.position3(time);
    double rho = std::sqrt(pos.X()*pos.X() + pos.Y()*pos.Y());
    double z = pos.Z();
    // radial velocity, signed by the extrapolation direction (positive = heading outward)
    double rvel = (rho > 0.0 ? (pos.X()*vel.X() + pos.Y()*vel.Y())/rho : 0.0)*timeDirSign(tdir);
    double vz = vel.Z()*timeDirSign(tdir); // z velocity, signed by the extrapolation direction
    if(debug_ > 2)std::cout << "R extrap start time " << time << " rho " << rho << " rvel " << rvel << " z " << z
      << " vz " << vz << " rval " << rval_ << " z[" << minz_ << "," << maxz_ << "]" << std::endl;
    // already outside the z bounds, or already at/beyond the target radius: stop
    if(z > maxz_ || z < minz_) return false;
    if(rho >= rval_) return false;
    // heading inward: this direction will never reach the outer target radius
    if(rvel <= 0.0) return false;
    // Decide GEOMETRICALLY which boundary the track reaches first (the BField domain step out here is set
    // by the field tolerance, not MaxDtStep, and can be large -- so a predicate that only stops AFTER the
    // step overshoots the z bound by a whole domain). A near-axial cosmic (large |v_z|, small v_rho)
    // reaches the z end before the target radius; extending it toward the radius is exactly what overshoots
    // the z end, carrying rho far past the DS shells and stranding them as interior pieces for the material
    // pass. In that case stop at the tracker boundary now and let the DS-material pass (no z bound, marches
    // shell-by-shell) extend it. Otherwise the track exits radially: keep going to the target radius.
    double time_to_rval = (rval_ - rho)/rvel;
    double time_to_zbound = std::numeric_limits<double>::max();
    if(vz > 0.0)      time_to_zbound = (maxz_ - z)/vz;
    else if(vz < 0.0) time_to_zbound = (z - minz_)/(-vz);
    if(time_to_zbound < time_to_rval) return false; // exits via the z end: stop at the tracker boundary
    return true; // radial exit: keep extending toward the target radius
  }
}
#endif
