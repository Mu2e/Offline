//
// Instantiate CLHEP::Hep3Vectors using polar coordinates.
// (Not provided in the the native class).
//
// $Id: polar3Vector.cc,v 1.3 2010/05/18 21:16:38 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/05/18 21:16:38 $
//
// Original author Rob Kutschke
//

#include <cmath>

#include "Mu2eUtilities/inc/ThreeVectorUtil.hh"
#include "Mu2eUtilities/inc/safeSqrt.hh"
#include "Mu2eUtilities/inc/sqrtOrThrow.hh"

using CLHEP::Hep3Vector;

namespace mu2e {

  CLHEP::Hep3Vector polar3Vector( double p0,
                                  double cz,
                                  double phi){

    double sz = safeSqrt( 1.-cz*cz);
    return CLHEP::Hep3Vector( p0*sz*cos(phi),
                              p0*sz*sin(phi),
                              p0*cz
                              );
  }

  CLHEP::Hep3Vector polar3Vector( double p0,
                                  double cz,
                                  double phi,
                                  double eps
                                  ){

    double sz = sqrtOrThrow( 1.-cz*cz,eps);
    return CLHEP::Hep3Vector( p0*sz*cos(phi),
                              p0*sz*sin(phi),
                              p0*cz
                              );
  }
  
} // end namespace mu2e
