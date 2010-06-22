//
// Manage all of the magnetic field maps for Mu2e.
//
// $Id: BFieldManager.cc,v 1.1 2010/06/22 16:44:25 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/06/22 16:44:25 $
//

// Includes from C++
#include <iostream>

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// Includes from Mu2e
#include "BFieldGeom/inc/BFieldManager.hh"

using namespace std;

namespace mu2e {

  BFieldManager::BFieldManager():
    _dsUniformValue(){
  }

  BFieldManager::~BFieldManager(){
  }

  // Get field at an arbitrary point.  Not yet implemented.
  // When written, this code will figure out which map to use
  // and will then look up the field in that map.
  CLHEP::Hep3Vector BFieldManager::getField( const CLHEP::Hep3Vector& point ) const{
    throw cms::Exception("GEOM")
      << "The getField method of BFieldManager is not yet implemented.\n"
      << "You can get a specific map and use the getField method of that map.\n";
  }

  // Get an arbitrary map, throw if it cannot be found.
  const BFMap& BFieldManager::getMapByName( const std::string key ) const{
    MapType::const_iterator i = _map.find(key);
    if ( i == _map.end() ){
      throw cms::Exception("GEOM")
        << "Requested map not found: " << key << "\n";      
    }
    return i->second;
  }

  // Create a new BFMap in the container of BFMaps.
  BFMap& BFieldManager::addBFMap( const std::string& key,
                                  CLHEP::Hep3Vector const& origin,
                                  int const nx, 
                                  int const ny, 
                                  int const nz ){

    // Construct an empty BFMap.
    BFMap bfmap(key, origin, nx, ny, nz);

    // Add it to the container of BFMaps.
    pair<MapType::iterator,bool> retval = _map.insert( MapType::value_type(key,bfmap) );

    // If there already was another Map with the same key, then it is a hard error.
    if ( !retval.second ){
      throw cms::Exception("GEOM")
        << "Trying to add a new magnetic field when the named field map already exists: "
        << key
        << "\n";      
    }

    // All Ok; return reference to the new BFMap.
    return retval.first->second;
  }


} // end namespace mu2e
