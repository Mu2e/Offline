//
// Manage all of the magnetic field maps for Mu2e.
//
// $Id: BFieldManager.cc,v 1.17 2013/08/30 22:25:58 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/08/30 22:25:58 $
//

// Includes from C++
#include <iostream>

// Framework includes
#include "cetlib_except/exception.h"

// Includes from Mu2e
#include "BFieldGeom/inc/BFieldManager.hh"

using namespace std;

namespace mu2e {

  // Get field at an arbitrary point. This code figures out which map to use
  // and looks up the field in that map.
  bool BFieldManager::getBFieldWithStatus(const CLHEP::Hep3Vector& point,
                                          CLHEP::Hep3Vector& result) const {

    const BFMap *m = cm_.findMap(point);

    if(m) {
      m->getBFieldWithStatus(point, result);
    }
    else {
      result = CLHEP::Hep3Vector(0.,0.,0.);
    }

    return (m != 0);
}

  bool BFieldManager::getNeighborPointBF (const CLHEP::Hep3Vector & testpoint, CLHEP::Hep3Vector neighborPoints[3], CLHEP::Hep3Vector neighborBF[3][3][3]) const {
    
    const BFMap *m = cm_.findMap(testpoint);

    if(m) {
      m->getNeighborPointBF(testpoint, neighborPoints, neighborBF);
    }
    else {
      for (int i = 0 ; i <3 ; ++i) {
        for (int j = 0 ; j <3 ; ++j) {
          for (int k = 0 ; k <3 ; ++k) {
            neighborBF[i][j][k] = CLHEP::Hep3Vector(0.,0.,0.);
          }
        }
        neighborPoints[i] = CLHEP::Hep3Vector(0., 0., 0.);
      }
    }

    return (m != 0);
  }



  // Create a new BFMap in the container of BFMaps.
  BFMap& BFieldManager::addBFMap(MapContainerType *mapContainer,
                                 const std::string& key,
                                 int nx, double xmin, double dx,
                                 int ny, double ymin, double dy,
                                 int nz, double zmin, double dz,
                                 BFMapType::enum_type type,
                                 double scaleFactor,
                                 BFInterpolationStyle interpStyle)
  {
    // If there already was another Map with the same key, then it is a hard error.
    if(!mapKeys_.insert(key).second) {
      throw cet::exception("GEOM")
        << "Trying to add a new magnetic field when the named field map already exists: "
        << key
        << "\n";
    }

    // Add an empty BFMap.
    mapContainer->push_back(BFMap(key, nx, xmin, dx, ny, ymin, dy, nz, zmin, dz, type, scaleFactor, interpStyle));

    return mapContainer->back();
  }

  void BFieldManager::print( ostream& out){
    out<<"================ BFieldManager: innerMaps ================\n";
    for ( MapContainerType::iterator i =innerMaps_.begin(); i != innerMaps_.end(); ++i ){
      i->print(out);
    }

    out<<"================ BFieldManager: outerMaps ================\n";
    for ( MapContainerType::iterator i =innerMaps_.begin(); i != innerMaps_.end(); ++i ){
      i->print(out);
    }

    out<<"================     BFieldManager end    ================\n";
  }

} // end namespace mu2e
