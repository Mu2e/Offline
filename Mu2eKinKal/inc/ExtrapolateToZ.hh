// predicate for extrapolation to a given Z position in a given direction
// track changes direction.
#ifndef Mu2eKinKal_ExtrapolateToZ_hh
#define Mu2eKinKal_ExtrapolateToZ_hh
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/General/TimeDir.hh"
#include "KinKal/Geometry/Plane.hh"
#include <limits>
#include "cetlib_except/exception.h"
namespace mu2e {
  using KinKal::TimeDir;
  using KinKal::timeDirSign;
  class ExtrapolateToZ {
    public:
      ExtrapolateToZ(double maxdt, double maxdtstep, double tol, double intertol, double zval,int debug) : maxDt_(maxdt), maxDtStep_(maxdtstep),
      dptol_(tol), intertol_(intertol), zval_(zval), debug_(debug), plane_(VEC3(0.0,0.0,1.0), VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,zval_)) {}
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
      double intertol_ = 1e-6; // intersection distance tolerance
      double zval_ = 0; // z value targeted
      int debug_ = 0; // debug level
      KinKal::Plane plane_; // KinKal plane
  };

  template <class KTRAJ> bool ExtrapolateToZ::needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, KinKal::TimeDir tdir) const {
    auto const& ktraj = tdir == TimeDir::forwards ? fittraj.back() : fittraj.front();
    auto time = tdir == TimeDir::forwards ? ktraj.range().end() : ktraj.range().begin();
    auto vel = ktraj.velocity(time)*timeDirSign(tdir); // sign by extrapolation direction
    auto epos =ktraj.position3(time);
    double dz = zval_ - epos.Z();
    if(debug_ > 2)std::cout << "Z extrap start time " << time << " dz " << dz << " zvel " << vel.Z() << " dz " <<  dz << std::endl;
    // false if we're heading away from the target z, otherwise true
    bool retval (true);
    if(dz*vel.Z() < 0) {
      // near reflection vz becomes an unreliable estimator of direction. Use the linearized direction instead
      auto axis = ktraj.linearize(time);
      double dist;
      plane_.intersect(axis,dist,tdir==TimeDir::forwards,intertol_);
      retval = dist*timeDirSign(tdir) > 0.0;
      if(debug_ > 2)std::cout << "Z extrap ray intersection result " << retval << " distance " << dist << std::endl;
    } else if (debug_ > 2)
      std::cout << "Z heading towards surface, distance " << dz << std::endl;
    return retval;
  }
}
#endif
