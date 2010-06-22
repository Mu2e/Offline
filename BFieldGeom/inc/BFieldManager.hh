#ifndef BFieldManager_HH
#define BFieldManager_HH

//
// Manage all of the magnetic field maps for Mu2e.
//
// $Id: BFieldManager.hh,v 1.1 2010/06/22 16:44:25 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/06/22 16:44:25 $
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
#include <map>

// Includes from Mu2e
#include "GeometryService/inc/Detector.hh"
#include "BFieldGeom/inc/BFMap.hh"

namespace mu2e {

  class BFieldManager : public Detector{

  public:

    // The class that knows how to populate this one.
    friend class BFieldManagerMaker;

    BFieldManager();
    ~BFieldManager();

    // Copying disabled - see below.

    virtual std::string name() const { return "BFieldManager";}

    // Get field at an arbitrary point.  Not yet implemented.
    CLHEP::Hep3Vector getField( const CLHEP::Hep3Vector& point ) const;

    // Get an arbitrary map.  Throw if it cannot be found.
    const BFMap& getMapByName( const std::string key ) const;

    // The uniform field in the DS is a special case.
    const CLHEP::Hep3Vector getDSUniformValue() const{
      return _dsUniformValue;
    }

    // Does this manager have a map with the given key.
    bool hasMapByName( const std::string key ) const{
      return ( _map.find(key) != _map.end() );
    }

  private:

    // Maps for various parts of the detector.
    typedef  std::map<std::string,BFMap> MapType;
    MapType _map;

    // Special case: uniform field in the DS.
    CLHEP::Hep3Vector _dsUniformValue;

    // Add an empty map to the list.  Used by BFieldManagerMaker.
    BFMap& addBFMap( const std::string& key,
                     CLHEP::Hep3Vector const& origin,
                     int const nx, 
                     int const ny, 
                     int const nz );

    // This class could support copying but it is not really needed and
    // I would like to prevent unintended copies ( people forgetting to
    // receive it into a const reference and making a copy instead ).
    BFieldManager( const BFieldManager& );
    BFieldManager& operator=( const BFieldManager& );

  }; // end class BFieldManager

} // end namespace mu2e
#endif
