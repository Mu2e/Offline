#ifndef ThreeVectorUtil_HH
#define ThreeVectorUtil_HH
//
// Various utilities that construct or operate on Hep3Vectors;
//
// $Id: ThreeVectorUtil.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//
// Notes
// 1) Hep3Vector itself provides setters for the representations:
//    ( r, theta, phi), ( r, phi, z), (r, eta, phi).
//    It does not provide ( r, cos(theta), phi).
//
// 2) Hep3Vector does not have a virtual destructor so we will
//    not inherit from it and will pay the copy penalty.
//
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e { 
  // Construct a Hep3Vector using the represenation:
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

