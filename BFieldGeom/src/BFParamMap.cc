//
// Class to hold one magnetic field map. The map
// is defined as a parametric function.
//
// $Id: BFParamMap.cc,v 1.00 2017/11/10 bpollack Exp $
// $Author: bpollack $
// $Date: 2017/11/10 $
//
// Original Brian Pollack, based on work by Krzysztof Genser, Rob Kutschke, Julie Managan, Bob
// Bernstein.

// C++ includes
#include <iomanip>
#include <iostream>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "BFieldGeom/inc/BFParamMap.hh"

// Other includes
#include "CLHEP/Vector/ThreeVector.h"
#include "cetlib_except/exception.h"
#include "math.h"

using namespace std;

namespace mu2e {

    // function to determine if the point is in the map; take into account Y-symmetry
    bool BFParamMap::isValid(CLHEP::Hep3Vector const& point) const {
        if (point.x() < _xmin || point.x() > _xmax) {
            return false;
        }
        if (point.y() < _ymin || point.y() > _ymax) {
            return false;
        }
        if (point.z() < _zmin || point.z() > _zmax) {
            return false;
        }
        return true;
    }

    bool BFParamMap::getBFieldWithStatus(const CLHEP::Hep3Vector& testpoint,
                                         CLHEP::Hep3Vector& result) const {
        bool retval(false);

        retval = fitFunction(testpoint, result);

        /*
        } else {
            throw cet::exception("GEOM")
                << "Unrecognized option for interpolation into the BField: " << _interpStyle
                << "\n";
        }
        */
        result *= _scaleFactor;
        return retval;
    }

    bool BFParamMap::fitFunction(const CLHEP::Hep3Vector& p, CLHEP::Hep3Vector& result) const {
        // Check validity.  Return a zero field and optionally print a warning.
        if (p.x() > (800 - 3896) or p.x() < (-800 - 3896) or p.y() > 800 or p.y() < -800 or
            p.z() < 4500 or p.z() > 13500) {
            return false;
        }

        if (!isValid(p)) {
            if (_warnIfOutside) {
                mf::LogWarning("GEOM")
                    << "Point is outside of the valid region of the map: " << _key << "\n"
                    << "Point in input coordinates: " << p << "\n";
            }
            result = CLHEP::Hep3Vector(0, 0, 0);
            return false;
        }

        vector<double> b_vec = _fitFunc->mag_field_function(p.x() + 3896, p.y(), p.z(), true);
        result = CLHEP::Hep3Vector(b_vec[0], b_vec[1], b_vec[2]);
        return true;
    }

    void BFParamMap::print(std::ostream& os) const {
        os << "Magnetic Field Info: " << _key << endl;

        cout << "Range X:    " << _xmin << " : " << _xmax << "  Middle: " << (_xmin + _xmax) / 2.
             << endl;
        cout << "Range Y:    " << _ymin << " : " << _ymax << "  Middle: " << (_ymin + _ymax) / 2.
             << endl;
        cout << "Range Z:    " << _zmin << " : " << _zmax << "  Middle: " << (_zmin + _zmax) / 2.
             << endl;

        if (_warnIfOutside) {
            cout << "Will warn if outside of the valid region." << endl;
        } else {
            cout << "Will not warn if outside of the valid region." << endl;
        }
    }

}  // end namespace mu2e
