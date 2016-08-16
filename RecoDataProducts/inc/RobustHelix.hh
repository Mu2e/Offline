//
//  Parameters describing a low-momentum helix that has non-zero transverse momentum
//  Original Author Dave Brown (LBNL) 15/8/2016
//
#ifndef RecoDataProducts_RobustHelix_HH
#define RecoDataProducts_RobustHelix_HH
// CLHEP
#include "CLHEP/Vector/ThreeVector.h"
// Root
#include "Rtypes.h"
namespece mu2e {
  class HelixVal;
  class RobustHelix {
    public:
      RobustHelix(Hep3Vector const& center, Double_t radius, Double_t lambda, Double_t fz0);
      ~RobustHelix(){}
      // accessors
      Hep3Vector const& center() const { return _center; }
      Double_t radius() const { return _radius; }
      Double_t lambda() const { return _lambda; }
      Double_t fz0() const { return _fz0; }
      double helicity() const { return copysign(1.0,lambda); } // check sign FIXME!
      // convert to BaBar rep
      void convert(HelixVal & hval) const;
    private:
      CLHEP::Hep3Vector _center; // circle center.  The z coordinate is not a fit parameter (nominally 0) (mm)
      Double_t _radius;  // transverse radius of the helix (mm)
      Double_t _lambda; // dz/dphi (mm/radian)
      Double_t _fz0; // azimuth (phi) at the center z position (radians)
  };
}
#endif
