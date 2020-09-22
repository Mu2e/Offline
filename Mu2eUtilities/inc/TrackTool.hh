#ifndef Mu2eUtilities_TrackTool_hh
#define Mu2eUtilities_TrackTool_hh
//
//
//
// Original author Rob Kutschke
//

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class TrackTool{

  public:
    TrackTool( int pdgId,
               double  q,
               CLHEP::Hep3Vector const& pos,
               CLHEP::Hep3Vector const& mom,
               double                   bz,
               CLHEP::Hep3Vector const& ref
               );

    // Accept compiler supplied d'tor, copy c'tor and assignment operator.

    // Accessors
    double cu() const { return cu_;}
    double da() const { return da_;}
    double fi() const { return fi_;}
    double ct() const { return ct_;}
    double z0() const { return z0_;}

    double rho() const { return radCurv_; }
    double xc() const { return xc_; }
    double yc() const { return yc_; }

    double d0x() const { return d0x_; }
    double d0y() const { return d0y_; }

    double u0() const { return u0_;}
    double v0() const { return v0_;}

    double cx() const { return cx_; }
    double cy() const { return cy_; }
    double cz() const { return cz_; }

    CLHEP::Hep3Vector center() const { return CLHEP::Hep3Vector(xc_, yc_, z0_); }

    // Turning angle from PCA to a given z.
    double psiAtZ( double z) const;

    // Position and momentum at a given turning angle.
    CLHEP::Hep3Vector positionAtPsi( double psi) const;
    CLHEP::Hep3Vector momentumAtPsi( double psi) const;

    // Position and momentum at a given z.
    CLHEP::Hep3Vector positionAtZ( double z) const;
    CLHEP::Hep3Vector momentumAtZ( double z) const;

  private:

    // PDG particle id.
    int pdgId_;

    // Electric charge of the track, in units of the proton charge.
    double q_;

    // A point on the track and the 3-momentum at that point:
    CLHEP::Hep3Vector pos_, p_;

    // Magnetic component field in the z direction, in T.
    double bz_;

    // Reference point.
    CLHEP::Hep3Vector ref_;

    double qq_;

    double cu_, da_, fi_, ct_, z0_;

    double s0_;

    double radCurv_;
    double pt_;
    double xc_, yc_;

    // Unit vector along track at xy plane
    double u0_, v0_;

    // Unit vector along track in 3D.
    double cx_, cy_, cz_;

    // Debug
    double d0x_, d0y_;

  };

} // namespace mu2e

#endif /* Mu2eUtilities_TrackTool_hh */
