#include "ExtinctionMonitorFNAL/Reconstruction/inc/TrackExtrapolator.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALMagnet.hh"

#include "RecoDataProducts/inc/ExtMonFNALTrkParam.hh"

#include <cassert>
#include <iostream>

namespace mu2e {
  namespace ExtMonFNAL {

    TrackExtrapolator::TrackExtrapolator(const ExtMon *detectorGeom)
      : extmon_(detectorGeom)
    {}

    //================================================================
    void TrackExtrapolator::extrapolateStraightLine(double newz0, ExtMonFNALTrkParam *pp) const {
      const double newx = pp->posx() + pp->slopex() * (newz0 - pp->z0());
      const double newy = pp->posy() + pp->slopey() * (newz0 - pp->z0());
      pp->setposx(newx);
      pp->setposy(newy);
      pp->setz0(newz0);
    }

    //================================================================
    void TrackExtrapolator::extrapolateToMagnet(ExtMonFNALTrkParam *pp) const {
      // FIXME:
      assert(pp->z0() > 0);

      const ExtMonFNALMagnet& mag = extmon_->spectrometerMagnet();
      const double nominalHalfBendAngle = mag.trackBendHalfAngle(mag.nominalMomentum());
      const double tana = tan(nominalHalfBendAngle);

      // distance from the ref point to stack axis crossing the magnet entrance plane
      const double z0entrance = mag.outerHalfSize()[2] / cos(nominalHalfBendAngle);

      // z coordinate at which the straight track enters the magnet
      const double zentrance = z0entrance + tana *
        ( (pp->z0() - z0entrance)*pp->slopey() - pp->posy()) / (1. + tana * pp->slopey());

      extrapolateStraightLine(zentrance, pp);
    }

    //================================================================
    // Input pars in the upstream stack system
    // output in the downstream stack system
    bool TrackExtrapolator::extrapolateThroughMagnet(ExtMonFNALTrkParam *pp) const {
      // FIXME: we assume an upstream-to-downstream extrapolation
      assert(pp->z0() > 0);

      const ExtMonFNALMagnet& mag = extmon_->spectrometerMagnet();
      const double nominalHalfBendAngle = mag.trackBendHalfAngle(mag.nominalMomentum());
      const double L = 2.*mag.outerHalfSize()[2];

      double sinSystem(0), cosSystem(0);
      sincos(nominalHalfBendAngle, &sinSystem, &cosSystem);
      //const double cosSystem = cos(nominalHalfBendAngle);
      //const double cosSystem = cos(nominalHalfBendAngle);

      // Compute track entrance coordinates in the magnet system: origin
      // at the ref point, Z parallel to the axis of the aperture, B
      // field along X.  That is, the system is rotated w.r.t. the UP
      // stack system by nominalHalfBendAngle in the YZ plane.
      //
      // Note that zEntrance==L/2, zExit==-L/2 where L is the magnet length.

      const double yEntrance = - sinSystem * pp->z0() + cosSystem * pp->posy();

      const double betaStack = atan(pp->slopey());
      const double beta = betaStack - nominalHalfBendAngle;

      if(beta <= -M_PI/2.) { // too steep angle, the track would go in the wrong direction
        return false;
      }

      const double tanbeta = tan(beta);
      const double cosbeta = cos(beta);

      const double xi = L * pp->rinv()/cosbeta;
      const double kappa = xi + 2*tanbeta;

      const double sqrtarg = 1- xi*kappa;
      if(sqrtarg < 0.) { // track does not reach the back of the magnet
        return false;
      }

      const double sqrtxikappa = sqrt(1-xi*kappa);

      // In the magnet system
      const double zExit = -L/2.;
      const double yExit = yEntrance - L*kappa/(1+sqrtxikappa);

      // In the downstream stack system: the same rotation as for UP to magnet
      const double zOutStack =   cosSystem * zExit + sinSystem * yExit;
      const double yOutStack = - sinSystem * zExit + cosSystem * yExit;

      // chord length in the projection on YZ plane
      const double chordLength = L* sqrt(1 + std::pow(kappa/(1+sqrtxikappa),2));

      // half bend angle for this track
      const double gamma = asin(xi*cosbeta*chordLength/(2*L));

      // The exit angle, in the magnet system
      const double betaExit = beta + 2*gamma;

      // arcLength = chordLength * gamma/sin(gamma)
      //
      // Using Euler's reflection formula
      // [Abramowitz & Stegun 6.1.17
      // after simple manipulations we get
      //
      const double arcLength = chordLength * tgamma(1.+gamma/M_PI) * tgamma(1.-gamma/M_PI);

      const double xOutStack = pp->posx() - arcLength*pp->slopex();
      const double syOutStack = tan(betaExit - nominalHalfBendAngle);

      // sx and rinv do not change
      pp->setz0(zOutStack);
      pp->setposx(xOutStack);
      pp->setposy(yOutStack);
      pp->setslopey(syOutStack);

      return true;
    }

    //================================================================
    bool TrackExtrapolator::extrapolateToPlane(unsigned plane, ExtMonFNALTrkParam *pp) const {
      bool res = true;
      const bool originDN(pp->z0()<0);
      const bool tgtDN(plane < extmon_->dn().nplanes());
      if(originDN ^ tgtDN) { // need to switch between stacks
        extrapolateToMagnet(pp);
        res = extrapolateThroughMagnet(pp);
        if(res) {
          const ExtMonFNALPlaneStack& stack = tgtDN ? extmon_->dn() : extmon_->up();
          const unsigned stackPlane = tgtDN ? plane : plane - extmon_->dn().nplanes();
          extrapolateStraightLine(stack.plane_zoffset()[stackPlane], pp);
        }
      }
      else { // stay in the same stack - just do straight line extrapolation
        const ExtMonFNALPlaneStack& stack = tgtDN ? extmon_->dn() : extmon_->up();
        const unsigned stackPlane = tgtDN ? plane : plane - extmon_->dn().nplanes();
        extrapolateStraightLine(stack.plane_zoffset()[stackPlane], pp);
      }
      return res;
    }

  } // namespace ExtMonFNAL
} // namespace mu2e
