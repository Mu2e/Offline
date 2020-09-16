//
// Manage all of the magnetic field maps for Mu2e.
//
//
// Modified by Brian Pollack to allow for polymorphic BField class.

// Includes from C++
#include <iostream>

// Framework includes
#include "cetlib_except/exception.h"

// Includes from Mu2e
#include "BFieldGeom/inc/BFieldManager.hh"

using namespace std;

namespace mu2e {

    // Get field at an arbitrary point. This code figures out which map to use
    // and looks up the field in that map.
    bool BFieldManager::getBFieldWithStatus(const CLHEP::Hep3Vector& point,
                                            CLHEP::Hep3Vector& result) const {
        return getBFieldWithStatus(point, cm_, result);
    }


    // Get field at an arbitrary point. This code figures out which map to use
    // and looks up the field in that map.
    bool BFieldManager::getBFieldWithStatus(const CLHEP::Hep3Vector& point,
                                            BFCacheManager const& cmgr,
                                            CLHEP::Hep3Vector& result) const {
        auto m = cmgr.findMap(point);

        if (m) {
            m->getBFieldWithStatus(point, result);
        } else {
            result = CLHEP::Hep3Vector(0., 0., 0.);
        }

        return (m != 0);
    }


    std::shared_ptr<BFGridMap> BFieldManager::addBFGridMap(MapContainerType* mapContainer,
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
                                                           BFInterpolationStyle interpStyle) {
        // If there already was another Map with the same key, then it is a hard error.
        if (!mapKeys_.insert(key).second) {
            throw cet::exception("GEOM")
                << "Trying to add a new magnetic field when the named field map already exists: "
                << key << "\n";
        }

        // Add an empty BFMap.
        auto new_map = std::make_shared<BFGridMap>(key, nx, xmin, dx, ny, ymin, dy, nz, zmin, dz,
                                                   type, scaleFactor, interpStyle);
        mapContainer->push_back(new_map);

        return new_map;
    }

    // Create a new BFGridMap in the container of BFMaps.
    std::shared_ptr<BFParamMap> BFieldManager::addBFParamMap(MapContainerType* mapContainer,
                                                             const std::string& key,
                                                             double xmin,
                                                             double xmax,
                                                             double ymin,
                                                             double ymax,
                                                             double zmin,
                                                             double zmax,
                                                             BFMapType::enum_type type,
                                                             double scaleFactor) {
        // If there already was another Map with the same key, then it is a hard error.
        if (!mapKeys_.insert(key).second) {
            throw cet::exception("GEOM")
                << "Trying to add a new magnetic field when the named field map already exists: "
                << key << "\n";
        }

        // Add an empty BFMap.
        auto new_map = std::make_shared<BFParamMap>(key, xmin, xmax, ymin, ymax, zmin, zmax, type,
                                                    scaleFactor);
        mapContainer->push_back(new_map);

        return new_map;
    }

    void BFieldManager::print(ostream& out) {
        out << "================ BFieldManager: innerMaps ================\n";
        for (MapContainerType::iterator i = innerMaps_.begin(); i != innerMaps_.end(); ++i) {
            (*i)->print(out);
        }

        out << "================ BFieldManager: outerMaps ================\n";
        for (MapContainerType::iterator i = outerMaps_.begin(); i != outerMaps_.end(); ++i) {
            (*i)->print(out);
        }

        out << "================     BFieldManager end    ================\n";
    }

}  // end namespace mu2e
