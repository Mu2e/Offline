//
//  $Id: TrkExtProtonAbsorber.cc,v 1.2 2012/08/18 22:27:56 mjlee Exp $
//  $Author: mjlee $
//  $Date: 2012/08/18 22:27:56 $
//
//  Original author MyeongJae Lee
//
//

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.

#include "CLHEP/Vector/ThreeVector.h"
#include "TrkExt/inc/TrkExtProtonAbsorber.hh"
#include "MECOStyleProtonAbsorberGeom/inc/MECOStyleProtonAbsorber.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeneralUtilities/inc/safeSqrt.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

using namespace CLHEP;

using namespace std;

namespace mu2e {




 
  TrkExtProtonAbsorber::TrkExtProtonAbsorber() :
    TrkExtShape(0.001),
    TrkExtMaterial("PE")
  {  
    valid = false;
  }

  void TrkExtProtonAbsorber::initialize() 
  {
    // geometry is read from mu2e coordinate. therefore, it should be transformed to detector coordinate
    art::ServiceHandle<GeometryService> geom;
    SimpleConfig const * config = &(geom->config());
    if (!(config->getBool("hasProtonAbsorber"))) {
      cerr << "TrkExtProtonAbsorber is disabled!" << endl;
      valid = false;
      return;
    }

    GeomHandle<MECOStyleProtonAbsorber> pabs;
    if (!(pabs->isAvailable(0)) && !(pabs->isAvailable(0)) ) {
      cerr << "TrkExtProtonAbsorber warning : pabs1 & pabs2 not avaliable" << endl;
      valid = false;
      return;
    }

    int index;
    if (pabs->isAvailable(0))  index = 0; 
    else    index = 1;
    z0    = pabs->part(index).center().z() - pabs->part(index).halfLength();
    r0in  = pabs->part(index).innerRadiusAtStart();
    r0out = pabs->part(index).outerRadiusAtStart();

    if (pabs->isAvailable(1)) index = 1; 
    else    index = 0;
    z1    = pabs->part(index).center().z() + pabs->part(index).halfLength();
    r1in  = pabs->part(index).innerRadiusAtEnd();
    r1out = pabs->part(index).outerRadiusAtEnd();

    GeomHandle<DetectorSystem> det;
    Hep3Vector origin = det->toMu2e( CLHEP::Hep3Vector(0.,0.,0.) );
    z0 -= origin.z();
    z1 -= origin.z();
    slope_out = (r1out - r0out) / (z1 - z0);
    slope_in  = (r1in - r0in) / (z1 - z0);

    cout << "TrkExtProtonAbsorber read : " << endl;
    cout << "  rOut = [" << r0out << ", " << r1out << "]" << endl;
    cout << "  rIn  = [" << r0in << ", " << r1in << "]" <<  endl;
    cout << "  z    = [" << z0 << ", " << z1 << "]" << endl;

    valid = true;
  
  }

  bool TrkExtProtonAbsorber::contains (Hep3Vector &xx) {
    // xx should be detector coordinate
    if (!valid) return false;
    double r = safeSqrt(xx.x() * xx.x() + xx.y() * xx.y());
    double z = xx.z();
    if (z <z0 || z1 < z) return false;
    if (r < slope_in*(z-z0)+r0in) return false;
    if (r > slope_out*(z-z0)+r0out) return false;
    return true;
  }





} // end namespace mu2e

