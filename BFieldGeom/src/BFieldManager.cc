//
// Manage all of the magnetic field maps for Mu2e.
//
// $Id: BFieldManager.cc,v 1.6 2011/03/08 00:40:46 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2011/03/08 00:40:46 $
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


  // Get field at an arbitrary point. This code figures out which map to use
  // and looks up the field in that map.
  bool BFieldManager::getBFieldWithStatus( const CLHEP::Hep3Vector & point,
					   CLHEP::Hep3Vector & result) const{

    // First check cached map
    if( _last_map!=0 && _last_map->getBFieldWithStatus(point,result) ) return true;

    // Loop over all maps and try to calculate field. 
    for( MapType::const_iterator im = _map.begin(); im!=_map.end(); ++im ) {
      if( im->second.getBFieldWithStatus(point,result) ) {
	_last_map = &(im->second);
	return true;
      }
    }
    
    _last_map = 0;
    result = CLHEP::Hep3Vector(0.,0.,0.);
    return false;

  }

  // Get an arbitrary map, throw if it cannot be found.
  // Search includes the manager class itself.
  const BFMapBase& BFieldManager::getMapByName( const std::string& key ) const{

    if( key == _key ) return (*this);
    return getContainedMapByName(key);

  }

  // Get an arbitrary map, throw if it cannot be found.
  // Search excludes the manager class itself.
  const BFMap& BFieldManager::getContainedMapByName( const std::string& key ) const{

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
                                  int const nz,
                                  BFMapType::enum_type type ){

    // Construct an empty BFMap.
    BFMap bfmap(key, nx, ny, nz, type);

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
  
  void BFieldManager::print( ostream& out){
    if ( _key.empty() ) {
      cout << "BFieldManager key name is the empty string." 
           << endl;
    } else{
      cout << "BFieldManager key name is the empty string." 
           << _key
           << endl;
    }
    for ( MapType::iterator i =_map.begin();
          i != _map.end(); ++i ){
      i->second.print(out);
    }
  }

} // end namespace mu2e
