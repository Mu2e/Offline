//
//  Parameters describing a low-momentum helix that has non-zero transverse momentum
//  Original Author Dave Brown (LBNL) 15/8/2016
//
// Mu2e
#include "RecoDataProducts/inc/RobustHelix.hh"
#include "RecoDataProducts/inc/HelixVal.hh"
#include "GeneralUtilities/inc/Angles.hh"
// C++
#include <math.h>
namespace mu2e {

  RobustHelix::RobustHelix(Hep3Vector const& center, Double_t radius, Double_t lambda, Double_t fz0) :
    _center(center), _radius(radius), _lambda(lambda), _fz0(fz0) {}

  bool RobustHelix::convert(HelixVal& hval) const {
    bool retval = _radius > 0.0;
// omega is the inverse transverse radius of the particle's circular motion.  Its
// signed by the helicity
//
    if(retval){

    // must specify the direction of time, that's implicit in the BaBar parameterization FIXME!!
      double omsign = copysign(1.0,-mytrk.particle().charge()*bz());
      hval._omega = omsign/radius;
      // phi0 is the azimuthal angle of the particle velocity vector at the point
      // of closest approach to the origin.  It's sign also depends on the angular
      // momentum.  To translate from the center, we need to reverse coordinates
      hval._phi0 = atan2(-omsign*_center.x(),omsign*_center.y());
      // d0 describes the distance to the origin at closest approach.
      // It is signed by the particle angular momentum WRT the origin.
      // The Helix fit radial bias is anti-correlated with d0; correct for it here.
      hval._d0 = omsign*(_center.perp() - _radius);
      // the dip angle is measured WRT the perpendicular
      hval._tanDip = omsign*_lambda/radius;
      // must change conventions here: fz0 is the phi at z=0, z0 is defined at the point of closest approach
      // resolve the loop ambiguity such that the POCA is closest to z=0.
      double dphi = Angles::deltaPhi(_fz0+omsign*M_PI_2,hval._phi0);
      // choose z0 (which loop) so that f=0 is as close to z=0 as possible
      hval._z0 = dphi*hval._tanDip/hval._omega;
    }
    return retval;
  }
}
