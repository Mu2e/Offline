#ifndef BFieldGeom_BFieldManager_hh
#define BFieldGeom_BFieldManager_hh
//
// Manage all of the magnetic field maps for Mu2e.
// This class holds the actual field maps, and provides an interface to compute B field.
//
// $Id: BFieldManager.hh,v 1.22 2013/08/30 22:25:22 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/08/30 22:25:22 $
//
// Notes:
// 1) This is a "dumb data" class. It does not know how to construct itself.
// 2) BFieldManagerMaker is the class that can populate this one.
// 3) Geant4 should access field maps via this class.
//

// C++ includes
#include <string>
#include <set>

// Includes from Mu2e
#include "Mu2eInterfaces/inc/Detector.hh"
#include "BFieldGeom/inc/BFMapType.hh"
#include "BFieldGeom/inc/BFMap.hh"
#include "BFieldGeom/inc/BFCacheManager.hh"
#include "BFieldGeom/inc/BFInterpolationStyle.hh"

namespace mu2e {

  class BFieldManagerMaker;

  class BFieldManager : virtual public Detector {

  public:

    // The class that knows how to populate this one.
    friend class BFieldManagerMaker;

    // Maps for various parts of the detector.
    typedef  std::vector<BFMap> MapContainerType;

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

    const MapContainerType& getInnerMaps() const { return innerMaps_; }
    MapContainerType& getInnerMaps() { return innerMaps_; }

    const MapContainerType& getOuterMaps() const { return outerMaps_; }
    MapContainerType& getOuterMaps() { return outerMaps_; }

    void print( std::ostream& out );

    bool getNeighborPointBF (const CLHEP::Hep3Vector &, CLHEP::Hep3Vector neighborPoints[3], CLHEP::Hep3Vector neighborBF[3][3][3]) const;


  private:
    // Private ctr.  An instance of BFieldManager should be obtained
    // via the BFieldManagerMaker class.
    BFieldManager() {}

    // This class could support copying but it is not really needed and
    // I would like to prevent unintended copies ( people forgetting to
    // receive it into a const reference and making a copy instead ).
    BFieldManager( const BFieldManager& );
    BFieldManager& operator=( const BFieldManager& );

    // Make sure map names are unique on the union of inner and outer maps
    std::set<std::string> mapKeys_;

    MapContainerType innerMaps_;
    MapContainerType outerMaps_;

    // Add an empty map to the list.  Used by BFieldManagerMaker.
    BFMap& addBFMap(MapContainerType *whichMap,
                    const std::string& key,
                    int nx, double xmin, double dx,
                    int ny, double ymin, double dy,
                    int nz, double zmin, double dz,
                    BFMapType::enum_type type,
                    double scaleFactor,
                    BFInterpolationStyle interpStyle );

    // Handles caching and overlap resolution logic
    BFCacheManager cm_;

  }; // end class BFieldManager

} // end namespace mu2e
#endif /* BFieldGeom_BFieldManager_hh */
