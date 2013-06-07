//
//  $Id: TrkExtToyDS.cc,v 1.5 2013/06/07 17:43:30 knoepfel Exp $
//  $Author: knoepfel $
//  $Date: 2013/06/07 17:43:30 $
//
//  Original author MyeongJae Lee
//
//

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.

#include "CLHEP/Vector/ThreeVector.h"
#include "TrkExt/inc/TrkExtToyDS.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "StoppingTargetGeom/inc/TargetFoil.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeneralUtilities/inc/safeSqrt.hh"

using namespace CLHEP;

using namespace std;

namespace mu2e {




 
  TrkExtToyDS::TrkExtToyDS() :
    TrkExtShape(0.001),
    TrkExtMaterial("Vac")
  { 
    rin = 1000;
    rout = 1300;
    zmin = 6000;
    zmax = 9226;
  }

  void TrkExtToyDS::initialize() {
    // geometry is read from mu2e coordinate. therefore, it should be transformed to detector coordinate
    GeomHandle<DetectorSystem> det;
    Hep3Vector origin = det->toMu2e( CLHEP::Hep3Vector(0.,0.,0.) );

    GeomHandle<DetectorSolenoid> dsgeom;
    rin = dsgeom->rIn1();
    rout = dsgeom->rOut2();
    zmin = dsgeom->position().z() - dsgeom->halfLength() - origin.z();
    zmax = dsgeom->position().z() + dsgeom->halfLength() - origin.z();

    cout << "TrkExtToyDS read rin=" << rin << ", rout=" << rout << ", zmin=" << zmin << ", zmax=" << zmax <<endl;
  }

  bool TrkExtToyDS::contains (Hep3Vector &xx) {
    // xx should be detector coordinate
    double r = safeSqrt(xx.x() * xx.x() + xx.y() * xx.y());
    double z = xx.z();
    if (z<zmin || z>zmax || r>rin) return false;
    return true;
  }





} // end namespace mu2e

