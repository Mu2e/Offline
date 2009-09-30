//
// Return Hep3Vector objects that are unit vectors uniformly
// distributed over the unit sphere.
// 
// $Id: RandomUnitSphere.cc,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

#include "CLHEP/Random/RandFlat.h"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/ThreeVectorUtil.hh"

using CLHEP::Hep3Vector;
using CLHEP::RandFlat;

namespace mu2e{ 

  Hep3Vector RandomUnitSphere::shoot() const{
    double  cz = _czmin  + ( _czmax  - _czmin  )*RandFlat::shoot();
    double phi = _phimin + ( _phimax - _phimin )*RandFlat::shoot();
    return polar3Vector ( 1., cz, phi);
  }

}
