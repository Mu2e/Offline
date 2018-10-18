//
// Return CLHEP::Hep3Vector objects that are unit vectors uniformly
// distributed over the unit sphere.
//
// $Id: RandomUnitSphere.cc,v 1.9 2012/02/03 05:02:43 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/03 05:02:43 $
//
// Original author Rob Kutschke
//

#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/ThreeVectorUtil.hh"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

using CLHEP::Hep3Vector;
using CLHEP::RandFlat;

namespace mu2e{

  RandomUnitSphere::RandomUnitSphere( CLHEP::HepRandomEngine& engine,
                                      double czmin,
                                      double czmax,
                                      double phimin,
                                      double phimax):
    _czmin(czmin),
    _czmax(czmax),
    _phimin(phimin),
    _phimax(phimax),
    _randFlat( engine ){
  }

  RandomUnitSphere::RandomUnitSphere( CLHEP::HepRandomEngine& engine,
                                      const RandomUnitSphereParams& pars):
    _czmin(pars.czmin),
    _czmax(pars.czmax),
    _phimin(pars.phimin),
    _phimax(pars.phimax),
    _randFlat( engine ){
  }

  CLHEP::Hep3Vector RandomUnitSphere::fire(){
    double  cz = _czmin  + ( _czmax  - _czmin  )*_randFlat.fire();
    double phi = _phimin + ( _phimax - _phimin )*_randFlat.fire();
    return polar3Vector ( 1., cz, phi);
  }

}
