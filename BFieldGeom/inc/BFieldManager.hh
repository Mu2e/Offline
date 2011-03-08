#ifndef BFieldManager_HH
#define BFieldManager_HH
//
// Manage all of the magnetic field maps for Mu2e.
//
// $Id: BFieldManager.hh,v 1.6 2011/03/08 00:40:23 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2011/03/08 00:40:23 $
//
// Notes:
// 1) This is a "dumb data" class. It does not know how to construct itself.
// 2) BFieldManagerMaker is the class that can populate this one.
// 3) Geant4 should access field maps via this class.
// 4) The case of a uniform BField in the DS is treated as a special case,
//    with special accessors.  Maybe it should not be treated this way but
//    it is for now.
//

// C++ includes
#include <string>
#include <iosfwd> 
#include <map>

// Includes from Mu2e
#include "GeometryService/inc/Detector.hh"
#include "BFieldGeom/inc/BFMapBase.hh"
#include "BFieldGeom/inc/BFMapType.hh"
#include "BFieldGeom/inc/BFMap.hh"

namespace mu2e {

  class BFieldManager : public Detector, public BFMapBase {

  public:

    // The class that knows how to populate this one.
    friend class BFieldManagerMaker;

    BFieldManager();
    ~BFieldManager();

    // Copying disabled - see below.

    std::string name() const { return "BFieldManager";}

    // Get field at an arbitrary point.
    bool getBFieldWithStatus(const CLHEP::Hep3Vector &, CLHEP::Hep3Vector& ) const;

    // Check if point belongs to any map
    bool isValid(CLHEP::Hep3Vector const& point) const;
    
    // The identifier for treating the BFieldManager as a field map itself.
    const std::string& getKey() const { return _key; };

    // Get an arbitrary map.  Throw if it cannot be found.
    const BFMapBase& getMapByName( const std::string& key ) const;

    // Get one of the contained maps (not the manager itself).  Throw if it cannot be found.
    // Return the BFMap interface, not the BFMapBase interface.
    const BFMap& getContainedMapByName( const std::string& key ) const;

    // The uniform field in the DS is a special case.
    const CLHEP::Hep3Vector getDSUniformValue() const{
      return _dsUniformValue;
    }

    // Does this manager have a map with the given key.
    bool hasMapByName( const std::string key ) const{
      return ( _map.find(key) != _map.end() );
    }

    double  rTorus() const { return _rTorus; }
    double xOffset() const { return _xOffset; }
    BFMapType type() const { return _type; }

    void print( std::ostream& out );

  private:

    // Name of the manager
    std::string _key;

    // Maps for various parts of the detector.
    typedef  std::map<std::string,BFMap> MapType;
    MapType _map;

    // Special case: uniform field in the DS.
    CLHEP::Hep3Vector _dsUniformValue;

    // Geometric properties of the muon beamline:
    // bend radius of the TS arcs and the offset from the center to DS or PS.
    double _rTorus;
    double _xOffset;

    // GMC, G4BL or possible future types.
    BFMapType _type;

    // Add an empty map to the list.  Used by BFieldManagerMaker.
    BFMap& addBFMap( const std::string& key,
                     int const nx, 
                     int const ny, 
                     int const nz,
                     BFMapType::enum_type type );

    // This class could support copying but it is not really needed and
    // I would like to prevent unintended copies ( people forgetting to
    // receive it into a const reference and making a copy instead ).
    BFieldManager( const BFieldManager& );
    BFieldManager& operator=( const BFieldManager& );

    // Use this pointer to cache access to maps collection
    mutable const BFMap * _last_map;

  }; // end class BFieldManager

} // end namespace mu2e
#endif
