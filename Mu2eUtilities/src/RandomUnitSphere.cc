//
// Return CLHEP::Hep3Vector objects that are unit vectors uniformly
// distributed over the unit sphere.
// 
// $Id: RandomUnitSphere.cc,v 1.2 2010/05/17 21:47:32 genser Exp $
// $Author: genser $ 
// $Date: 2010/05/17 21:47:32 $
//
// Original author Rob Kutschke
//

#include "CLHEP/Random/RandFlat.h"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/ThreeVectorUtil.hh"

using CLHEP::Hep3Vector;
using CLHEP::RandFlat;

namespace mu2e{ 

  CLHEP::Hep3Vector RandomUnitSphere::shoot() const{
    double  cz = _czmin  + ( _czmax  - _czmin  )*CLHEP::RandFlat::shoot();
    double phi = _phimin + ( _phimax - _phimin )*CLHEP::RandFlat::shoot();
    return polar3Vector ( 1., cz, phi);
  }

}
