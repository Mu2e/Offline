//
// Class to hold one magnetic field map. The map
// is defined as a parametric function.
//
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
#include <gsl/gsl_sf_bessel.h>
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

        retval = evalFit(testpoint, result);

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

    bool BFParamMap::evalFit(const CLHEP::Hep3Vector& p, CLHEP::Hep3Vector& result) const {
        // Check validity.  Return a zero field and optionally print a warning.
        if (!isValid(p)) {
            if (_warnIfOutside) {
                mf::LogWarning("GEOM")
                    << "Point is outside of the valid region of the map: " << _key << "\n"
                    << "Point in input coordinates: " << p << "\n";
            }
            result = CLHEP::Hep3Vector(0, 0, 0);
            return false;
        }

        double r, phi;
        // The declarations below are to reduce computation time
        double bessels[2];
        double tmp_rho;
        double cos_nphi, cos_kmsz;
        double sin_nphi, sin_kmsz;
        double abp, abm;
        vector<double> out(3, 0);
        phi = atan2(p.y(), p.x() + 3896);
        r = sqrt(pow(p.x() + 3896, 2) + pow(p.y(), 2));
        double abs_r = abs(r);

        for (int n = 0; n < _ns; ++n) {
            for (int m = 1; m <= _ms; ++m) {
                tmp_rho = _kms[n][m - 1] * abs_r;
                bessels[0] = gsl_sf_bessel_In(n, tmp_rho);
                bessels[1] = gsl_sf_bessel_In(n + 1, tmp_rho);
                _iv[n][m - 1] = bessels[0];
                if (tmp_rho == 0) {
                    _ivp[n][m - 1] = 0.5 * (gsl_sf_bessel_In(n - 1, 0) + bessels[1]);
                } else {
                    _ivp[n][m - 1] = (n / tmp_rho) * bessels[0] + bessels[1];
                }
            }
        }

        double br(0.0);
        double bphi(0.0);
        double bz(0.0);
        // Here is the meat of the calculation:
        for (int n = 0; n < _ns; ++n) {
            cos_nphi = cos(n * phi + _Ds[n]);
            sin_nphi = -sin(n * phi + _Ds[n]);
            for (int m = 0; m < _ms; ++m) {
                cos_kmsz = cos(_kms[n][m] * p.z());
                sin_kmsz = sin(_kms[n][m] * p.z());
                abp = _As[n][m] * cos_kmsz + _Bs[n][m] * sin_kmsz;
                abm = -_As[n][m] * sin_kmsz + _Bs[n][m] * cos_kmsz;
                br += cos_nphi * _ivp[n][m] * _kms[n][m] * abp;
                bz += cos_nphi * _iv[n][m] * _kms[n][m] * abm;
                if (abs_r > 1e-10) {
                    bphi += n * sin_nphi * (1 / abs_r) * _iv[n][m] * abp;
                }
            }
        }

        double cp = cos(phi);
        double sp = sin(phi);
        result = CLHEP::Hep3Vector(br * cp - bphi * sp, br * sp + bphi * cp, bz);

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

    void BFParamMap::calcConstants() {
        for (int n = 0; n < _ns; ++n) {
            _kms.push_back(vector<double>());
            for (int m = 1; m <= _ms; ++m) {
                _kms[n].push_back(m * M_PI / _Reff);
            }
        }

        for (int n = 0; n < _ns; ++n) {
            vector<double> tmp_v1(_ms, 0.0);
            vector<double> tmp_v2(_ms, 0.0);
            _iv.push_back(tmp_v1);
            _ivp.push_back(tmp_v2);
        }
    }

}  // end namespace mu2e
