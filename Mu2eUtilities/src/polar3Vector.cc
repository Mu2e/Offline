//
// Instantiate Hep3Vectors using polar coordinates.
// (Not provided in the the native class).
//
// $Id: polar3Vector.cc,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

#include <cmath>

#include "Mu2eUtilities/inc/ThreeVectorUtil.hh"
#include "Mu2eUtilities/inc/safeSqrt.hh"
#include "Mu2eUtilities/inc/sqrtOrThrow.hh"

using CLHEP::Hep3Vector;

namespace mu2e {

  Hep3Vector polar3Vector( double p0,
			   double cz,
			   double phi){

    double sz = safeSqrt( 1.-cz*cz);
    return Hep3Vector( p0*sz*cos(phi),
		       p0*sz*sin(phi),
		       p0*cz
		       );
  }

  Hep3Vector polar3Vector( double p0,
			   double cz,
			   double phi,
			   double eps
			   ){

    double sz = sqrtOrThrow( 1.-cz*cz,eps);
    return Hep3Vector( p0*sz*cos(phi),
		       p0*sz*sin(phi),
		       p0*cz
		       );
  }
  
}
