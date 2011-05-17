//
// Return CLHEP::Hep3Vector objects that are unit vectors uniformly
// distributed over the unit sphere.
// 
// $Id: RandomUnitSphere.cc,v 1.6 2011/05/17 15:36:01 greenc Exp $
// $Author: greenc $ 
// $Date: 2011/05/17 15:36:01 $
//
// Original author Rob Kutschke
//

#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/ThreeVectorUtil.hh"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/RandomNumberGeneratorService.h"

using CLHEP::Hep3Vector;
using CLHEP::RandFlat;

namespace mu2e{ 

  RandomUnitSphere::RandomUnitSphere( double czmin,
                                      double czmax,
                                      double phimin,
                                      double phimax):
    _czmin(czmin),
    _czmax(czmax),
    _phimin(phimin),
    _phimax(phimax),
    _randFlat( art::ServiceHandle<art::RandomNumberGeneratorService>()->getEngine() ){
  }

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

  CLHEP::Hep3Vector RandomUnitSphere::fire(){
    double  cz = _czmin  + ( _czmax  - _czmin  )*_randFlat.fire();
    double phi = _phimin + ( _phimax - _phimin )*_randFlat.fire();
    return polar3Vector ( 1., cz, phi);
  }

}
