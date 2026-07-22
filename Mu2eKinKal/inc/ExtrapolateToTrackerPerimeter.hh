// Predicate for extrapolation to the tracker PERIMETER: the closed cylinder bounding the tracker
// active volume, {rho = sqrt(x^2+y^2) <= R} intersected with {zmin <= z <= zmax}. The track is
// extended until it reaches the FIRST bounding face it meets -- the curved outer wall (constant rho)
// OR either z-end plane (constant z) -- whichever comes first for that track.
//
// This is the direction-agnostic generalization of the two former predicates: ExtrapolateToZ (a
// downstream track exits through a constant-z end plane) and ExtrapolateToR (a cosmic traverses the
// tracker radially -- mostly in y -- and exits through the curved outer wall). One predicate handles
// both, choosing the face per track from its own kinematics, so no configuration switch is needed.
//
// All three bounds (R, zmin, zmax) come from the tracker geometry (KinKalGeom::tracker(), itself
// sourced from GeometryService), so the tracker perimeter is defined in exactly one place -- no
// hand-entered radius/z limits that must be kept in sync with the geometry.
#ifndef Mu2eKinKal_ExtrapolateToTrackerPerimeter_hh
#define Mu2eKinKal_ExtrapolateToTrackerPerimeter_hh
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/General/TimeDir.hh"
#include <limits>
#include <cmath>
#include <iostream>
#include "cetlib_except/exception.h"
namespace mu2e {
  using KinKal::TimeDir;
  using KinKal::timeDirSign;
  class ExtrapolateToTrackerPerimeter {
    public:
      ExtrapolateToTrackerPerimeter(double maxdt, double maxdtstep, double tol,
          double rval, double zmin, double zmax, int debug) :
        maxDt_(maxdt), maxDtStep_(maxdtstep), dptol_(tol), rval_(rval), zmin_(zmin), zmax_(zmax), debug_(debug) {}
      // interface for extrapolation
      double maxDt() const { return maxDt_; } // maximum time to extend the track, WRT the time of the first(last) measurement
      double maxDtStep() const { return maxDtStep_; }
      double dpTolerance() const { return dptol_; } // tolerance on fractional momentum change
      double rVal() const { return rval_; }
      double zMin() const { return zmin_; }
      double zMax() const { return zmax_; }
      int debug() const { return debug_; }
      // extrapolation predicate: the track will be extrapolated until this predicate returns false, subject to the maximum time
      template <class KTRAJ> bool needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const;
    private:
      double maxDt_ = -1; // maximum extrapolation time
      double maxDtStep_ = -1; // maximum extrapolation time step in a single iteration
      double dptol_ = 1e10; // fractional momentum tolerance in BField domain
      double rval_ = 0; // outer wall transverse radius (mm)
      double zmin_ = -std::numeric_limits<double>::max(); // upstream z-end plane (mm)
      double zmax_ =  std::numeric_limits<double>::max(); // downstream z-end plane (mm)
      int debug_ = 0; // debug level
  };

  template <class KTRAJ> bool ExtrapolateToTrackerPerimeter::needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, KinKal::TimeDir tdir) const {
    auto const& ktraj = tdir == TimeDir::forwards ? fittraj.back() : fittraj.front();
    auto time = tdir == TimeDir::forwards ? ktraj.range().end() : ktraj.range().begin();
    auto vel = ktraj.velocity(time);
    auto pos = ktraj.position3(time);
    double rho = std::sqrt(pos.X()*pos.X() + pos.Y()*pos.Y());
    double z = pos.Z();
    // velocities signed by the extrapolation direction (positive = heading toward that boundary)
    double rvel = (rho > 0.0 ? (pos.X()*vel.X() + pos.Y()*vel.Y())/rho : 0.0)*timeDirSign(tdir); // outward radial
    double vz = vel.Z()*timeDirSign(tdir);
    if(debug_ > 2)std::cout << "Perimeter extrap start time " << time << " rho " << rho << " rvel " << rvel
      << " z " << z << " vz " << vz << " R " << rval_ << " z[" << zmin_ << "," << zmax_ << "]" << std::endl;
    // already at/outside any bounding face: stop (the recording step finds the precise crossing)
    if(rho >= rval_) return false;
    if(z <= zmin_ || z >= zmax_) return false;
    // time to reach each face we are actually approaching (INF if heading away from it)
    static const double INF = std::numeric_limits<double>::max();
    double t_r  = rvel > 0.0 ? (rval_ - rho)/rvel : INF; // outer wall
    double t_zlo = vz < 0.0  ? (z - zmin_)/(-vz)  : INF; // upstream z-end
    double t_zhi = vz > 0.0  ? (zmax_ - z)/vz     : INF; // downstream z-end
    // inside the volume but heading toward no face (e.g. a stalled step): nothing more to reach
    if(t_r == INF && t_zlo == INF && t_zhi == INF) return false;
    // still strictly inside and closing on a face: keep extending
    return true;
  }
}
#endif
