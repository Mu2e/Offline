//
//
//  Original author MyeongJae Lee
//
//

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.

#include "CLHEP/Vector/ThreeVector.h"
#include "Offline/TrkExt/inc/TrkExtToyDS.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/StoppingTargetGeom/inc/StoppingTarget.hh"
#include "Offline/StoppingTargetGeom/inc/TargetFoil.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/GeneralUtilities/inc/safeSqrt.hh"

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

