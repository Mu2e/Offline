//
// G4 interface to the Detector Solenoid full magnetic field.
//
//
// Original author Julie Managan and Bob Bernstein
// Major rewrite by Rob Kutschke at version 1.4

// C++ includes
//#include <cmath>
#include <iostream>

// Mu2e includes.
#include "Mu2eG4/inc/Mu2eGlobalField.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "BFieldGeom/inc/BFieldManager.hh"

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/ThreeVector.h"


#include "G4Threading.hh"

using CLHEP::Hep3Vector;
using namespace std;

namespace mu2e {

  Mu2eGlobalField::Mu2eGlobalField(const G4ThreeVector& mapOrigin)
  {
    // Load map.
    update(mapOrigin);
  }

  // This is the entry point called by G4.
  void Mu2eGlobalField::GetFieldValue(const G4double Point[4],
                              G4double *Bfield) const {

    // Put point in required format and required reference frame.
    CLHEP::Hep3Vector point(Point[0],Point[1],Point[2]);
    point -= _mapOrigin;

    // Look up BField and reformat to required return format.
    const CLHEP::Hep3Vector bf = _map->getBField(point, _cm);
    Bfield[0] = bf.x()*CLHEP::tesla;
    Bfield[1] = bf.y()*CLHEP::tesla;
    Bfield[2] = bf.z()*CLHEP::tesla;

    /*
    cout << "Mu2eGlobalField map=" << _map->getKey()
         << " point=("<<point.x()<<","<<point.y()<<","<<point.z()<<")"
         << " field=("<<Bfield[0]<<","<<Bfield[1]<<","<<Bfield[2]<<")"
         << endl;
    */

  }

  // Update the map and its origin.  Might be called for new runs?
  void Mu2eGlobalField::update( const G4ThreeVector& mapOrigin){

    _mapOrigin = mapOrigin;

    // Handle to the BField manager.
    GeomHandle<BFieldManager> bfMgr;

    // Throws if the map is not found.
    _map = &*bfMgr;

    _cm = bfMgr->cacheManager();
      
      //std::cout << " from Thread #" << G4Threading::G4GetThreadId()
      //<< ", address of CacheManager is " << &_cm << std::endl;
  }

} // end namespace mu2e
