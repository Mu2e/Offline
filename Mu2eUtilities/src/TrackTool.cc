// Convert track from 3-point + 3-momentum representation to other
//representations.
//
//
// Original author Rob Kutschke
//

#include <cmath>
#include <iostream>

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/ThreeVector.h"


#include "Mu2eUtilities/inc/TrackTool.hh"

using namespace std;
//using CLHEP::Hep3Vector;

namespace mu2e {

  TrackTool::TrackTool( int pdgId,
                        double                   q,
                        CLHEP::Hep3Vector const& pos,
                        CLHEP::Hep3Vector const& p,
                        double                   bz,
                        CLHEP::Hep3Vector const& ref
                        ):
    pdgId_(pdgId),
    q_(q),
    pos_(pos),
    p_(p),
    bz_(bz),
    ref_(ref){

    // The conversion factor from pT to curvature.
    // Bz is signed so a carries "geometric charge".
    static double c_b = CLHEP::c_light*1.e-3;
    double          a = c_b * bz_ * q_;

    // The geoemetric charge.
    qq_ = ( a > 0. ) ? 1. : -1.;

    // Position relative to the reference point.
    double dx=pos.x()-ref.x();
    double dy=pos.y()-ref.y();
    double dz=pos.z()-ref.z();

    // Quantities that depend only on the momentum components.
    pt_      = p.perp();
    cu_      = 0.5 * a / pt_;
    ct_      = p.z() / pt_;
    radCurv_ = std::abs(pt_/a);
    cx_      = p.x()/p.mag();
    cy_      = p.y()/p.mag();
    cz_      = p.z()/p.mag();

    // Center of curvature and unsigned radius from center of curvature to the reference point
    xc_ = dx + qq_*p.y()/pt_*radCurv_;
    yc_ = dy - qq_*p.x()/pt_*radCurv_;
    double rc  = sqrt( xc_*xc_ + yc_*yc_);

    // Unit vector to center of curvature from reference point.
    double xhc = xc_/rc;
    double yhc = yc_/rc;

    // Point of closest approach, in 2D, to the reference point.
    d0x_ = xc_ - xhc*radCurv_;
    d0y_ = yc_ - yhc*radCurv_;

    // Azimuth, its cosine and sine.
    // Need to check if this works for positive geometric charge.
    double f = ( radCurv_ > rc )? -qq_ : qq_;
    fi_ = atan2( f*d0x_, -f*d0y_ );
    u0_ = cos(fi_);
    v0_ = sin(fi_);

    // Signed impact parameter.
    da_ = -d0x_*v0_ + d0y_*u0_;

    // psi is the 2D turning angle from the PCA to the input point;
    // psi lies on [0,2pi).
    // cpsi and spsi are its cosine and sine
    double cpsi =  qq_*( xhc*p.y()/pt_ - yhc*p.x()/pt_);
    double spsi = -qq_*( xhc*p.x()/pt_ + yhc*p.y()/pt_);
    double psi  = atan2(spsi,cpsi);
    if ( psi < 0. ) psi += 2.*M_PI;

    // Z0 parameter at PCA.
    z0_ = dz - radCurv_*ct_*psi;

    /*
    // Check the result.
    //CLHEP::Hep3Vector qmom(momentumAtPsi(psi));
    CLHEP::Hep3Vector qmom(momentumAtZ(pos.z()));

    cout << "px: "
         << qmom.x()    << " "
         << p.x() << " "
         << qmom.x()-p.x()
         << endl;

    cout << "py: "
         << qmom.y()    << " "
         << p.y() << " "
         << qmom.y()-p.y()
         << endl;

    cout << "pz: "
         << qmom.z()    << " "
         << p.z() << " "
         << qmom.z()-p.z()
         << endl;

    //CLHEP::Hep3Vector xyz(positionAtPsi(psi));
    CLHEP::Hep3Vector xyz(positionAtZ(pos.z()));
    cout << "x: "
         << xyz.x() << " "
         << pos.x() << " "
         << xyz.x()-pos.x()
         << endl;

    cout << "y: "
         << xyz.y() << " "
         << pos.y() << " "
         << xyz.y()-pos.y()
         << endl;

    cout << "z: "
         << xyz.z() << " "
         << pos.z() << " "
         << xyz.z()-pos.z()
         << endl;

    */
  } // end TrackTool::TrackTool

  // Turning angle from PCA to a given z.
  double TrackTool::psiAtZ( double z) const{
    return (z - z0_ )/radCurv_/ct_;
  }

  CLHEP::Hep3Vector TrackTool::positionAtPsi( double psi) const{
    double x = xc_ - qq_*radCurv_*sin( fi_ - qq_*psi);
    double y = yc_ + qq_*radCurv_*cos( fi_ - qq_*psi);
    double z = z0_ + radCurv_*psi*ct_;
    return CLHEP::Hep3Vector(x,y,z);
  }

  CLHEP::Hep3Vector TrackTool::momentumAtPsi( double psi) const{
    double px = pt_*cos(fi_ - qq_*psi);
    double py = pt_*sin(fi_ - qq_*psi);
    double pz = pt_*ct_;
    return CLHEP::Hep3Vector(px,py,pz);
  }

  CLHEP::Hep3Vector TrackTool::positionAtZ( double z) const{
    double psi = psiAtZ(z);
    return positionAtPsi(psi);
  }

  CLHEP::Hep3Vector TrackTool::momentumAtZ( double z) const{
    double psi = psiAtZ(z);
    return momentumAtPsi(psi);
  }


}
