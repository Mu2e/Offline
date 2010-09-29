//
// Manage all of the magnetic field maps for Mu2e.
//
// $Id: BFieldManager.cc,v 1.3 2010/09/29 22:51:43 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/09/29 22:51:43 $
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
    _key(""), _dsUniformValue(), _last_map(0) {
  }

  BFieldManager::~BFieldManager(){
  }

  // Check if point belongs to any map
  bool BFieldManager::isValid(CLHEP::Hep3Vector const& point) const {

    // First check cached map
    if( _last_map!=0 && _last_map->isValid(point) ) return true;

    // Loop over all maps
    for( MapType::const_iterator im = _map.begin(); im!=_map.end(); ++im ) {
      if( im->second.isValid(point) ) {
	_last_map = &(im->second);
	return true;
      }
    }

    return false;

  }


  // Get field at an arbitrary point.  Not yet implemented.
  // When written, this code will figure out which map to use
  // and will then look up the field in that map.
  CLHEP::Hep3Vector BFieldManager::getBField( const CLHEP::Hep3Vector& point ) const{

    // First check cached map
    if( _last_map!=0 && _last_map->isValid(point) ) return _last_map->getBField(point);

    // Loop over all maps and check if point belongs to them. Use first found map.
    for( MapType::const_iterator im = _map.begin(); im!=_map.end(); ++im ) {
      if( im->second.isValid(point) ) {
	_last_map = &(im->second);
	return im->second.getBField(point);
      }
    }

    return CLHEP::Hep3Vector(0.,0.,0.);

  }

  // Get an arbitrary map, throw if it cannot be found.
  const BFMapBase& BFieldManager::getMapByName( const std::string key ) const{

    if( key == _key ) return (*this);

    MapType::const_iterator i = _map.find(key);
    if ( i == _map.end() ){
      throw cms::Exception("GEOM")
        << "Requested map not found: " << key << "\n";      
    }
    return i->second;

  }

  // Create a new BFMap in the container of BFMaps.
  BFMap& BFieldManager::addBFMap( const std::string& key,
                                  int const nx, 
                                  int const ny, 
                                  int const nz ){

    // Construct an empty BFMap.
    BFMap bfmap(key, nx, ny, nz);

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
