#ifndef ThreeVectorUtil_HH
#define ThreeVectorUtil_HH
//
// Various utilities that construct or operate on CLHEP::Hep3Vectors;
//
// $Id: ThreeVectorUtil.hh,v 1.2 2010/05/17 21:47:32 genser Exp $
// $Author: genser $ 
// $Date: 2010/05/17 21:47:32 $
//
// Original author Rob Kutschke
//
// Notes
// 1) CLHEP::Hep3Vector itself provides setters for the representations:
//    ( r, theta, phi), ( r, phi, z), (r, eta, phi).
//    It does not provide ( r, cos(theta), phi).
//
// 2) CLHEP::Hep3Vector does not have a virtual destructor so we will
//    not inherit from it and will pay the copy penalty.
//
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e { 
  // Construct a CLHEP::Hep3Vector using the represenation:
  // ( magnitude, cos(polar angle), azimuth) = ( r, cz, phi)
  // Uses safesqrt to compute sin(theta).
  CLHEP::Hep3Vector polar3Vector( double r,
				  double cz,
				  double phi );

  // As above but uses a tolerance in the computation of sin(theta).
  // ( uses sqrtOrThrow ).
  CLHEP::Hep3Vector polar3Vector( double r,
				  double cz,
				  double phi,
				  double eps );
}

#endif

