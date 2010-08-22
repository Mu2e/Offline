//
// Return CLHEP::Hep3Vector objects that are unit vectors uniformly
// distributed over the unit sphere.
// 
// $Id: RandomUnitSphere.cc,v 1.4 2010/08/22 20:09:52 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/08/22 20:09:52 $
//
// Original author Rob Kutschke
//

#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/ThreeVectorUtil.hh"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/RandomNumberGeneratorService.h"

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
    _randFlat( edm::Service<edm::RandomNumberGeneratorService>()->getEngine() ){
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

  CLHEP::Hep3Vector RandomUnitSphere::fire() const{
    double  cz = _czmin  + ( _czmax  - _czmin  )*_randFlat.fire();
    double phi = _phimin + ( _phimax - _phimin )*_randFlat.fire();
    return polar3Vector ( 1., cz, phi);
  }

  CLHEP::Hep3Vector RandomUnitSphere::fire( double magnitude ) const{
    return magnitude*fire();
  }

}
