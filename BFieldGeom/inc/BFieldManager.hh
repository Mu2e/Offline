#ifndef BFieldGeom_BFieldManager_hh
#define BFieldGeom_BFieldManager_hh
//
// Manage all of the magnetic field maps for Mu2e.
// This class holds the actual field maps, and provides an interface to compute B field.
//
// Modified by Brian Pollack to allow for polymorphic BField class.
//
// Notes:
// 1) This is a "dumb data" class. It does not know how to construct itself.
// 2) BFieldManagerMaker is the class that can populate this one.
// 3) Geant4 should access field maps via this class.
//

// C++ includes
#include <set>
#include <string>

// Includes from Mu2e
#include "BFieldGeom/inc/BFCacheManager.hh"
#include "BFieldGeom/inc/BFGridMap.hh"
#include "BFieldGeom/inc/BFInterpolationStyle.hh"
#include "BFieldGeom/inc/BFMap.hh"
#include "BFieldGeom/inc/BFMapType.hh"
#include "BFieldGeom/inc/BFParamMap.hh"
#include "DataProducts/inc/XYZVec.hh"
#include "Mu2eInterfaces/inc/Detector.hh"

namespace mu2e {

    class BFieldManagerMaker;

    class BFieldManager : virtual public Detector {
       public:
        // The class that knows how to populate this one.
        friend class BFieldManagerMaker;

        // Maps for various parts of the detector.
        typedef std::vector<std::shared_ptr<BFMap>> MapContainerType;

        // Get field at an arbitrary point.
        bool getBFieldWithStatus(const CLHEP::Hep3Vector&, CLHEP::Hep3Vector&) const;
        bool getBFieldWithStatus(const CLHEP::Hep3Vector&,
                                 BFCacheManager const&,
                                 CLHEP::Hep3Vector&) const;

        // Just return zero for out of range.
        CLHEP::Hep3Vector getBField(const CLHEP::Hep3Vector& pos) const {
            // Default c'tor sets all components to zero - which is what we need here.
            CLHEP::Hep3Vector result;
            getBFieldWithStatus(pos, result);
            return result;
        }

        CLHEP::Hep3Vector getBField(const CLHEP::Hep3Vector& pos,
                                    BFCacheManager const& cmgr) const {
            // Default c'tor sets all components to zero - which is what we need here.
            CLHEP::Hep3Vector result;
            getBFieldWithStatus(pos, cmgr, result);
            return result;
        }

        XYZVec getBField(const XYZVec& pos) const {
          // Default c'tor sets all components to zero - which is what we need here.
          CLHEP::Hep3Vector b;

          CLHEP::Hep3Vector p(pos.x(), pos.y(), pos.z());
          getBFieldWithStatus(p, b);
          XYZVec result(b.x(), b.y(), b.z());

          return result;
        }

        BFCacheManager cacheManager() const { return cm_; }

        const MapContainerType& getInnerMaps() const { return innerMaps_; }
        MapContainerType& getInnerMaps() { return innerMaps_; }

        const MapContainerType& getOuterMaps() const { return outerMaps_; }
        MapContainerType& getOuterMaps() { return outerMaps_; }

        void print(std::ostream& out);

       private:
        // Private ctr.  An instance of BFieldManager should be obtained
        // via the BFieldManagerMaker class.
        BFieldManager() {}

        // This class could support copying but it is not really needed and
        // I would like to prevent unintended copies ( people forgetting to
        // receive it into a const reference and making a copy instead ).
        BFieldManager(const BFieldManager&);
        BFieldManager& operator=(const BFieldManager&);

        // Make sure map names are unique on the union of inner and outer maps
        std::set<std::string> mapKeys_;

        MapContainerType innerMaps_;
        MapContainerType outerMaps_;

        // Add an empty grid-like map to the list.  Used by BFieldManagerMaker.
        std::shared_ptr<BFGridMap> addBFGridMap(MapContainerType* whichMap,
                                                const std::string& key,
                                                int nx,
                                                double xmin,
                                                double dx,
                                                int ny,
                                                double ymin,
                                                double dy,
                                                int nz,
                                                double zmin,
                                                double dz,
                                                BFMapType::enum_type type,
                                                double scaleFactor,
                                                BFInterpolationStyle interpStyle);

        // Add an empty parametric map to the list.  Used by BFieldManagerMaker.
        std::shared_ptr<BFParamMap> addBFParamMap(MapContainerType* whichMap,
                                                  const std::string& key,
                                                  double xmin,
                                                  double xmax,
                                                  double ymin,
                                                  double ymax,
                                                  double zmin,
                                                  double zmax,
                                                  BFMapType::enum_type type,
                                                  double scaleFactor);

        // Handles caching and overlap resolution logic
        BFCacheManager cm_;

    };  // end class BFieldManager

}  // end namespace mu2e
#endif /* BFieldGeom_BFieldManager_hh */
