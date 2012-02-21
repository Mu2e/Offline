#ifndef BFieldGeom_BFieldManager_hh
#define BFieldGeom_BFieldManager_hh
//
// Manage all of the magnetic field maps for Mu2e.
//
// $Id: BFieldManager.hh,v 1.14 2012/02/21 22:26:23 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/21 22:26:23 $
//
// Notes:
// 1) This is a "dumb data" class. It does not know how to construct itself.
// 2) BFieldManagerMaker is the class that can populate this one.
// 3) Geant4 should access field maps via this class.
//

// C++ includes
#include <string>
#include <iosfwd>
#include <map>

// Includes from Mu2e
#include "GeometryService/inc/Detector.hh"
#include "BFieldGeom/inc/BFMapType.hh"
#include "BFieldGeom/inc/BFMap.hh"

namespace mu2e {

  class BFieldManagerMaker;

  class BFieldManager : public Detector {

  public:

    // The class that knows how to populate this one.
    friend class BFieldManagerMaker;

    // Maps for various parts of the detector.
    typedef  std::map<std::string,BFMap> MapContainerType;

    // Implement the Detector method
    virtual std::string name() const { return "BFieldManager";}

    // Get field at an arbitrary point.
    bool getBFieldWithStatus(const CLHEP::Hep3Vector &, CLHEP::Hep3Vector& ) const;

    // Just return zero for out of range.
    CLHEP::Hep3Vector getBField(const CLHEP::Hep3Vector& pos) const {
      CLHEP::Hep3Vector result;
      if(!getBFieldWithStatus(pos,result)) {
        result = CLHEP::Hep3Vector(0,0,0);
      }
      return result;  // make sure NRVO can be applied
    }

    // Check if point belongs to any map
    bool isValid(CLHEP::Hep3Vector const& point) const;

    // Get an arbitrary map.  Throw if it cannot be found.
    const BFMap& getMapByName( const std::string& key ) const;

    // Does this manager have a map with the given key.
    bool hasMapByName( const std::string key ) const{
      return ( _map.find(key) != _map.end() );
    }

    const MapContainerType& getMapContainer() const { return _map; }

    void print( std::ostream& out );

  private:
    // Private ctr.  An instance of BFieldManager should be obtained
    // via the BFieldManagerMaker class.
    BFieldManager();

    // This class could support copying but it is not really needed and
    // I would like to prevent unintended copies ( people forgetting to
    // receive it into a const reference and making a copy instead ).
    BFieldManager( const BFieldManager& );
    BFieldManager& operator=( const BFieldManager& );

    MapContainerType _map;

    // Add an empty map to the list.  Used by BFieldManagerMaker.
    BFMap& addBFMap( const std::string& key,
                     int const nx,
                     int const ny,
                     int const nz,
                     BFMapType::enum_type type,
                     double scaleFactor );

    // Use this pointer to cache access to maps collection
    mutable const BFMap * _last_map;

  }; // end class BFieldManager

} // end namespace mu2e
#endif /* BFieldGeom_BFieldManager_hh */
