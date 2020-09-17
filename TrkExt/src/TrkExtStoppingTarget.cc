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
#include "TrkExt/inc/TrkExtStoppingTarget.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "StoppingTargetGeom/inc/TargetFoil.hh"
#include "GeneralUtilities/inc/safeSqrt.hh"

using namespace CLHEP;

using namespace std;

namespace mu2e {




 
  TrkExtStoppingTarget::TrkExtStoppingTarget() :
    TrkExtShape(),
    TrkExtMaterial("Al")
  { 
   foil.clear();
   nfoil = -1;
  }

  void TrkExtStoppingTarget::initialize() {
    // geometry is read from mu2e coordinate. therefore, it should be transformed to detector coordinate
    GeomHandle<DetectorSystem> det;
    Hep3Vector origin = det->toMu2e( CLHEP::Hep3Vector(0.,0.,0.) );

    GeomHandle<StoppingTarget> target;
    foil.clear();
    for (int i = 0 ; i < target->nFoils() ; ++i) {
      TargetFoil f = target->foil(i);
      foil_data_type ftmp (f.rOut(),
                           f.centerInDetectorSystem().z() - f.halfThickness(),
                           f.centerInDetectorSystem().z(),
                           f.centerInDetectorSystem().z() + f.halfThickness());
      foil.push_back(ftmp);
    }

    nfoil = foil.size();

    rmax = -999;
    zmax = -999999;
    zmin = 999999;
    for (int i = 0 ; i < nfoil ; ++i) {
      cout << "TrkExtStoppingTarget read foil " << i << " at " << foil[i].zc << ", r = " << foil[i].rout << endl;
      if (foil[i].rout > rmax) rmax = foil[i].rout;
      if (foil[i].z1 > zmax) zmax = foil[i].z1;
      if (foil[i].z0 < zmin) zmin = foil[i].z0;
    }
    cout << "TrkExtStoppingTarget read foil rmax=" << rmax << ", zmax=" << zmax << ", zmin=" << zmin << endl;
  }

  bool TrkExtStoppingTarget::contains (Hep3Vector &xx) {
    // xx should be detector coordinate
    if (nfoil<=0) return false;
    double r = safeSqrt(xx.x() * xx.x() + xx.y() * xx.y());
    double z = xx.z();
    if (z<zmin || z>zmax || r>rmax) return false;
    for (int i = 0 ; i <nfoil ; ++i) {
      if (foil[i].z0 < z && z <foil[i].z1 && r <foil[i].rout) return true;
    }
    return false;
  }





} // end namespace mu2e

