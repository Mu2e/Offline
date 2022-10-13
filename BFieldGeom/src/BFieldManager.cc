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
#include "Offline/BFieldGeom/inc/BFieldManager.hh"

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


  BFieldManager::BFieldManager(MapContainerType const& innerMaps,
                               MapContainerType const& outerMaps):
    innerMaps_(innerMaps),outerMaps_(outerMaps) {
    cm_.setMaps(innerMaps, outerMaps);

  }
    void BFieldManager::print(ostream& out) const {
        out << "================ BFieldManager: innerMaps ================\n";
        for (auto i = innerMaps_.begin(); i != innerMaps_.end(); ++i) {
            (*i)->print(out);
        }

        out << "================ BFieldManager: outerMaps ================\n";
        for (auto i = outerMaps_.begin(); i != outerMaps_.end(); ++i) {
            (*i)->print(out);
        }

        out << "================     BFieldManager end    ================\n";
    }

}  // end namespace mu2e
