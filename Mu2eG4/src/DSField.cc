//
// G4 interface to the Detector Solenoid full magnetic field.
//
// $Id: DSField.cc,v 1.4 2010/06/22 16:42:22 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/06/22 16:42:22 $
//
// Original author Julie Managan and Bob Bernstein
// Major rewrite by Rob Kutschke at version 1.4

// C++ includes
//#include <cmath>
#include <iostream>

// Mu2e includes.
#include "Mu2eG4/inc/DSField.hh"
#include "BFieldGeom/inc/BFMap.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "BFieldGeom/inc/BFieldManager.hh"

// Other includes.
#include "CLHEP/Vector/ThreeVector.h"

using CLHEP::Hep3Vector;
using namespace std;

namespace mu2e {

  DSField::DSField( std::string name, G4ThreeVector mapOrigin ):
    _name(name){

    // Load map.
    update(mapOrigin);
  }

  // This is the entry point called by G4.
  void DSField::GetFieldValue(const G4double Point[4],
                              G4double *Bfield) const {

    // Put point in required format and required reference frame.
    CLHEP::Hep3Vector point(Point[0],Point[1],Point[2]);
    point -= _mapOrigin;

    // Look up BField and reformat to required return format.
    const CLHEP::Hep3Vector bf = _map->getBField(point);
    Bfield[0] = bf.x()*CLHEP::tesla;
    Bfield[1] = bf.y()*CLHEP::tesla;
    Bfield[2] = bf.z()*CLHEP::tesla;

  }

  // Update the map and its origin.  Might be called for new runs?
  void DSField::update( const G4ThreeVector& mapOrigin){

    _mapOrigin = mapOrigin;

    // Handle to the BField manager.
    GeomHandle<BFieldManager> bfMgr;

    // Throws if the map is not found.
    _map = &bfMgr->getMapByName(_name);
  }

} // end namespace mu2e
