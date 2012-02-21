#ifndef BFieldGeom_BFieldManager_hh
#define BFieldGeom_BFieldManager_hh
//
// Manage all of the magnetic field maps for Mu2e.
//
// $Id: BFieldManager.hh,v 1.13 2012/02/21 22:26:02 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/21 22:26:02 $
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
#include "BFieldGeom/inc/BFMapType.hh"
#include "BFieldGeom/inc/BFMap.hh"

namespace mu2e {

  class BFieldManager : public Detector {

  public:

    // The class that knows how to populate this one.
    friend class BFieldManagerMaker;

    // Maps for various parts of the detector.
    typedef  std::map<std::string,BFMap> MapContainerType;

    BFieldManager();
    ~BFieldManager();

    // Copying disabled - see below.

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

    // The uniform field in the DS is a special case.
    const CLHEP::Hep3Vector getDSUniformValue() const{
      return _dsUniformValue;
    }

    // The gradient field in the DS is a special case.
    const CLHEP::Hep3Vector getDSGradientValue() const{
      return _dsGradientValue;
    }

    // Does this manager have a map with the given key.
    bool hasMapByName( const std::string key ) const{
      return ( _map.find(key) != _map.end() );
    }

    double  rTorus() const { return _rTorus; }
    double xOffset() const { return _xOffset; }
    BFMapType type() const { return _type; }

    const MapContainerType& getMapContainer() const { return _map; }

    void print( std::ostream& out );

  private:

    MapContainerType _map;

    // Special case: uniform field in the DS.
    CLHEP::Hep3Vector _dsUniformValue;

    // Special case: gradient field in the DS
    CLHEP::Hep3Vector _dsGradientValue;

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
                     BFMapType::enum_type type,
                     double scaleFactor );

    // This class could support copying but it is not really needed and
    // I would like to prevent unintended copies ( people forgetting to
    // receive it into a const reference and making a copy instead ).
    BFieldManager( const BFieldManager& );
    BFieldManager& operator=( const BFieldManager& );

    // Use this pointer to cache access to maps collection
    mutable const BFMap * _last_map;

  }; // end class BFieldManager

} // end namespace mu2e
#endif /* BFieldGeom_BFieldManager_hh */
