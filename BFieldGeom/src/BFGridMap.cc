//
// Class to hold one magnetic field map. The map
// is defined on a regular cartesian grid.
//
//
// Original Rob Kutschke, based on work by Julie Managan and Bob Bernstein.
// Rewritten in part by Krzysztof Genser to correct mistake pointed by RB and to save execution
// time. Rewritten again by Brian Pollack to seperate out interpolation methods from parametric
// methods.

// C++ includes
#include <iomanip>
#include <iostream>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "Offline/BFieldGeom/inc/BFGridMap.hh"
#include "Offline/BFieldGeom/inc/BFMap.hh"

// Other includes
#include "CLHEP/Vector/ThreeVector.h"
#include "cetlib_except/exception.h"
#include <cmath>

using namespace std;

namespace mu2e {

    // function to determine if the point is in the map; take into account Y-symmetry
    bool BFGridMap::isValid(CLHEP::Hep3Vector const& point) const {
        if (point.x() < _xmin || point.x() > _xmax) {
            return false;
        }
        if (_flipy) {
            if (fabs(point.y()) < _ymin || fabs(point.y()) > _ymax) {
                return false;
            }
        } else {
            if (point.y() < _ymin || point.y() > _ymax) {
                return false;
            }
        }
        if (point.z() < _zmin || point.z() > _zmax) {
            return false;
        }

        return true;
    }

    CLHEP::Hep3Vector BFGridMap::cellFraction(const CLHEP::Hep3Vector& pos,
                                              const GridPoint& ipos) const {
        const CLHEP::Hep3Vector gridpos(grid2point(ipos.ix, ipos.iy, ipos.iz));
        return CLHEP::Hep3Vector((pos.x() - gridpos.x()) / _dx, (pos.y() - gridpos.y()) / _dy,
                                 (pos.z() - gridpos.z()) / _dz);
    }

    // Populate a 3x3x3 array with the field values centered on the grid point
    // that is the closest to the requested point.
    bool BFGridMap::getNeighbors(int ix,
                                 int iy,
                                 int iz,
                                 CLHEP::Hep3Vector neighborsBF[3][3][3]) const {
        for (int i = 0; i != 3; ++i) {
            unsigned int xindex = ix + i - 1;
            for (int j = 0; j != 3; ++j) {
                unsigned int yindex = iy + j - 1;
                for (int k = 0; k != 3; ++k) {
                    unsigned int zindex = iz + k - 1;
                    if (!_isDefined(xindex, yindex, zindex))
                        return false;
                    neighborsBF[i][j][k] = _field(xindex, yindex, zindex);
                    /*
                              cout << "Neighbor(" << xindex << "," << yindex << "," << zindex
                              << ") = (" << neighborsBF(i,j,k).x() << ","
                              << neighborsBF(i,j,k).y() << "," << neighborsBF(i,j,k).z()
                              << ")" << endl;
                    */
                }
            }
        }
        return true;
    }

    // Function to interpolate the BField value at the point from the values
    // at in the neighbor grid.
    CLHEP::Hep3Vector BFGridMap::interpolate(CLHEP::Hep3Vector const vec[3][3][3],
                                             const CLHEP::Hep3Vector& frac) const {
        // Create vecs and vectors
        double x1d[3], y1d[3], z1d[3];
        CLHEP::Hep3Vector vecx[9];
        CLHEP::Hep3Vector vecxy[3];
        double xin = frac[0];
        double yin = frac[1];
        double zin = frac[2];

        /*  cout << "--------Starting interpolation-------- \n" << "Input: ("
            << xin << "," << yin << "," << zin << ")" << endl;*/

        // First loop - interpolate xin
        for (int j = 0; j != 3; ++j) {
            for (int k = 0; k != 3; ++k) {
                for (int i = 0; i != 3; ++i) {
                    /*        cout << "1st loop access: (" << vec[i][j][k].x() << ","
                              << vec[i][j][k].y() << "," << vec[i][j][k].z() << ")" << endl;*/
                    x1d[i] = (vec[i][j][k]).x();
                    y1d[i] = (vec[i][j][k]).y();
                    z1d[i] = (vec[i][j][k]).z();
                }
                vecx[j * 3 + k] =
                    CLHEP::Hep3Vector(gmcpoly2(x1d, xin), gmcpoly2(y1d, xin), gmcpoly2(z1d, xin));
            }
        }

        // Second loop - interpolate yin
        for (int k = 0; k != 3; ++k) {
            for (int j = 0; j != 3; ++j) {
                int index = j * 3 + k;
                x1d[j] = vecx[index].x();
                y1d[j] = vecx[index].y();
                z1d[j] = vecx[index].z();
            }
            vecxy[k] =
                CLHEP::Hep3Vector(gmcpoly2(x1d, yin), gmcpoly2(y1d, yin), gmcpoly2(z1d, yin));
        }

        // Third loop - interpolate zin
        for (int k = 0; k != 3; ++k) {
            x1d[k] = vecxy[k].x();
            y1d[k] = vecxy[k].y();
            z1d[k] = vecxy[k].z();
        }

        // Return BField 3Vector
        return CLHEP::Hep3Vector(gmcpoly2(x1d, zin), gmcpoly2(y1d, zin), gmcpoly2(z1d, zin));
    }

    // Standard Lagrange formula for 2nd order polynomial fit of
    // univariate function
    double BFGridMap::gmcpoly2(double const f1d[3], double const& x) const {
        static const double x0(0.), x1(1.), x2(2.);

        return f1d[0] * (x - x1) * (x - x2) / ((x0 - x1) * (x0 - x2)) +
               f1d[1] * (x - x0) * (x - x2) / ((x1 - x0) * (x1 - x2)) +
               f1d[2] * (x - x0) * (x - x1) / ((x2 - x0) * (x2 - x1));
    }

    bool BFGridMap::getBFieldWithStatus(const CLHEP::Hep3Vector& testpoint,
                                        CLHEP::Hep3Vector& result) const {
        bool retval(false);

        if (_interpStyle == BFInterpolationStyle::trilinear) {
            retval = interpolateTriLinear(testpoint, result);

        } else {
            throw cet::exception("GEOM")
                << "Unrecognized option for interpolation into the BField: " << _interpStyle
                << "\n";
        }
        result *= _scaleFactor;
        return retval;
    }

    // The algorithm is:
    // Find the grid cube in which the point lives - this defines eight corner points.
    // Assign a weight to each corner that is the "distance" to each corner - see below for
    // its precise definition.  The field value at the test point is the weighted sum of
    // each of the 8 corner points.
    bool BFGridMap::interpolateTriLinear(const CLHEP::Hep3Vector& p,
                                         CLHEP::Hep3Vector& result) const {
        double px = p.x();
        double py = p.y();
        if (_flipy)
            py = std::abs(p.y());
        double pz = p.z();

        // Indicies into each dimension;
        int i = floor((px - _xmin) / _dx);
        int j = floor((py - _ymin) / _dy);
        int k = floor((pz - _zmin) / _dz);

        // Check that we are inside the map.
        if (i < 0 || i >= int(_nx) || j < 0 || j >= int(_ny) || k < 0 || k >= int(_nz)) {
            if (_warnIfOutside) {
                mf::LogWarning("GEOM")
                    << "Point is outside of the valid region of the map: " << _key << "\n"
                    << "Point in input coordinates: " << p << "\n";
            }
            result = CLHEP::Hep3Vector(0., 0., 0.);
            return false;
        }

        // Trilinear fractional weighting factors.
        double fx = 1.0 - (px - _xmin - i * _dx) / _dx;
        double fy = 1.0 - (py - _ymin - j * _dy) / _dy;
        double fz = 1.0 - (pz - _zmin - k * _dz) / _dz;

        // Field values at the 8 corner points.
        // Guess that a copy is faster than a pointer for reasons of locality
        // of reference in the downstream code?
        CLHEP::Hep3Vector c[8] = {_field(i, j, k),         _field(i + 1, j, k),
                                  _field(i, j + 1, k),     _field(i + 1, j + 1, k),
                                  _field(i, j, k + 1),     _field(i + 1, j, k + 1),
                                  _field(i, j + 1, k + 1), _field(i + 1, j + 1, k + 1)};

        double bx = c[0].x() * fx * fy * fz + c[1].x() * (1.0 - fx) * fy * fz +
                    c[2].x() * fx * (1.0 - fy) * fz + c[3].x() * (1.0 - fx) * (1.0 - fy) * fz +
                    c[4].x() * fx * fy * (1.0 - fz) + c[5].x() * (1.0 - fx) * fy * (1.0 - fz) +
                    c[6].x() * fx * (1.0 - fy) * (1.0 - fz) +
                    c[7].x() * (1.0 - fx) * (1.0 - fy) * (1.0 - fz);

        double by = c[0].y() * fx * fy * fz + c[1].y() * (1.0 - fx) * fy * fz +
                    c[2].y() * fx * (1.0 - fy) * fz + c[3].y() * (1.0 - fx) * (1.0 - fy) * fz +
                    c[4].y() * fx * fy * (1.0 - fz) + c[5].y() * (1.0 - fx) * fy * (1.0 - fz) +
                    c[6].y() * fx * (1.0 - fy) * (1.0 - fz) +
                    c[7].y() * (1.0 - fx) * (1.0 - fy) * (1.0 - fz);

        double bz = c[0].z() * fx * fy * fz + c[1].z() * (1.0 - fx) * fy * fz +
                    c[2].z() * fx * (1.0 - fy) * fz + c[3].z() * (1.0 - fx) * (1.0 - fy) * fz +
                    c[4].z() * fx * fy * (1.0 - fz) + c[5].z() * (1.0 - fx) * fy * (1.0 - fz) +
                    c[6].z() * fx * (1.0 - fy) * (1.0 - fz) +
                    c[7].z() * (1.0 - fx) * (1.0 - fy) * (1.0 - fz);

        // Need the signed value of p.y() here - the variable py will not do.
        if (_flipy && p.y() < 0)
            by = -by;

        result = CLHEP::Hep3Vector(bx, by, bz);

        return true;
    }


    bool BFGridMap::getNeighborPointBF(const CLHEP::Hep3Vector& testpoint,
                                       CLHEP::Hep3Vector neighborPoints[3],
                                       CLHEP::Hep3Vector neighborBF[3][3][3]) const {
        // Allow y-symmetry if grid is only defined for y > 0;
        int sign(1);
        CLHEP::Hep3Vector point(testpoint.x(), testpoint.y(), testpoint.z());
        if (_flipy && testpoint.y() < 0) {
            sign = -1;
            double y = -testpoint.y();
            point.setY(y);
        }

        // Check validity.  Return a zero field and optionally print a warning.
        if (!isValid(point)) {
            if (_warnIfOutside) {
                mf::LogWarning("GEOM")
                    << "Point is outside of the valid region of the map: " << _key << "\n"
                    << "Point in input coordinates: " << testpoint << "\n";
            }
            return false;
        }

        // Get the indices of the nearest grid point
        unsigned int ix = static_cast<int>((point.x() - _xmin) / _dx + 0.5);
        unsigned int iy = static_cast<int>((point.y() - _ymin) / _dy + 0.5);
        unsigned int iz = static_cast<int>((point.z() - _zmin) / _dz + 0.5);

        // Correct for edge points by moving their NGPt just inside the edge
        if (ix == 0) {
            ++ix;
        }
        if (ix == _nx - 1) {
            --ix;
        }
        if (iy == 0) {
            ++iy;
        }
        if (iy == _ny - 1) {
            --iy;
        }
        if (iz == 0) {
            ++iz;
        }
        if (iz == _nz - 1) {
            --iz;
        }

        // check if the point had a field defined

        if (!_isDefined(ix, iy, iz)) {
            if (_warnIfOutside) {
                mf::LogWarning("GEOM")
                    << "Point's field is not defined in the map: " << _key << "\n"
                    << "Point in input coordinates: " << testpoint << "\n";
                mf::LogWarning("GEOM") << "ix=" << ix << " iy=" << iy << " iz=" << iz << "\n";
            }
            return false;
        }

        // set neighborPoints and BFs

        for (int i = 0; i != 3; ++i) {
            unsigned int xindex = ix + i - 1;
            for (int j = 0; j != 3; ++j) {
                unsigned int yindex = iy + j - 1;
                for (int k = 0; k != 3; ++k) {
                    unsigned int zindex = iz + k - 1;
                    if (!_isDefined(xindex, yindex, zindex)) {
                        if (_warnIfOutside) {
                            mf::LogWarning("GEOM")
                                << "Point's neighboring field is not defined in the map: " << _key
                                << "\n"
                                << "Point in input coordinates: " << testpoint << "\n";
                            mf::LogWarning("GEOM")
                                << "ix=" << ix << " iy=" << iy << " iz=" << iz << "\n";
                        }
                        return false;
                    }
                    neighborBF[i][j][k] = _field(xindex, yindex, zindex);
                    // Reassign y sign
                    if (_flipy && sign == -1) {
                        neighborBF[i][j][k].setY(-neighborBF[i][j][k].y());
                    }
                }
            }
        }

        for (int i = 0; i < 3; ++i) {
            neighborPoints[i].set(double(ix + i - 1) * _dx + _xmin,
                                  double(iy + i - 1) * _dy + _ymin,
                                  double(iz + i - 1) * _dz + _zmin);
        }

        return true;
    }

    void BFGridMap::print(std::ostream& os) const {
        os << "Magnetic Field Info: " << _key << endl;

        cout << "Dimensions: " << _nx << " " << _ny << " " << _nz << endl;
        cout << "Range X:    " << _xmin << " : " << _xmax << "  Middle: " << (_xmin + _xmax) / 2.
             << endl;
        cout << "Range Y:    " << _ymin << " : " << _ymax << "  Middle: " << (_ymin + _ymax) / 2.
             << endl;
        cout << "Range Z:    " << _zmin << " : " << _zmax << "  Middle: " << (_zmin + _zmax) / 2.
             << endl;
        cout << "Distance:       " << _dx << " " << _dy << " " << _dz << endl;

        cout << "Field at the edges: " << _field(0, 0, 0) << ", " << _field(_nx - 1, 0, 0) << ", "
             << _field(0, _ny - 1, 0) << ", " << _field(0, 0, _nz - 1) << ", "
             << _field(_nx - 1, _ny - 1, 0) << ", " << _field(_nx - 1, _ny - 1, _nz - 1) << endl;

        cout << "Field in the middle: " << _field(_nx / 2, _ny / 2, _nz / 2) << endl;

        if (_warnIfOutside) {
            cout << "Will warn if outside of the valid region." << endl;
        } else {
            cout << "Will not warn if outside of the valid region." << endl;
        }
    }

    std::ostream& operator<<(std::ostream& os, const BFGridMap::GridPoint& p) {
        return os << "BFGridMap::GridPoint(" << p.ix << ", " << p.iy << ", " << p.iz << ")";
    }

}  // end namespace mu2e
