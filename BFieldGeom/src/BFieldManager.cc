//
// Manage all of the magnetic field maps for Mu2e.
//
// $Id: BFieldManager.cc,v 1.15 2012/02/29 00:34:28 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/29 00:34:28 $
//

// Includes from C++
#include <iostream>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

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

  // Create a new BFMap in the container of BFMaps.
  BFMap& BFieldManager::addBFMap(MapContainerType *mapContainer,
                                 const std::string& key,
                                 int nx, double xmin, double dx,
                                 int ny, double ymin, double dy,
                                 int nz, double zmin, double dz,
                                 BFMapType::enum_type type,
                                 double scaleFactor)
  {
    // If there already was another Map with the same key, then it is a hard error.
    if(!mapKeys_.insert(key).second) {
      throw cet::exception("GEOM")
        << "Trying to add a new magnetic field when the named field map already exists: "
        << key
        << "\n";
    }

    // Add an empty BFMap.
    mapContainer->push_back(BFMap(key, nx, xmin, dx, ny, ymin, dy, nz, zmin, dz, type, scaleFactor));

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
