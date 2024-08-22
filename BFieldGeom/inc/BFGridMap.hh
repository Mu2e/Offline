#ifndef BFieldGeom_BFGridMap_hh
#define BFieldGeom_BFGridMap_hh
//
// Class to hold one magnetic field map. The map is defined on a regular cartesian grid.
// All field maps are given in the standard Mu2e coordinate system.
// Units are: space point in mm, field values in tesla.
//
//
// Original Rob Kutschke, based on work by Julie Managan and Bob Bernstein.
// Rewritten in part by Krzysztof Genser to save execution time
// Rewritten again by Brian Pollack to separate out Grid-like maps from other map types
//

//#include <iosfwd>
#include <ostream>
#include <string>
#include "Offline/BFieldGeom/inc/BFInterpolationStyle.hh"
#include "Offline/BFieldGeom/inc/BFMap.hh"
#include "Offline/BFieldGeom/inc/BFMapType.hh"
#include "Offline/BFieldGeom/inc/Container3D.hh"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {
    class BFGridMap : public BFMap {
       public:
        friend class BFieldManagerMaker;

        struct GridPoint {
            unsigned ix;
            unsigned iy;
            unsigned iz;
            GridPoint(unsigned a, unsigned b, unsigned c) : ix(a), iy(b), iz(c) {}
        };

        BFGridMap(std::string const& filename,
                  int nx,
                  double xmin,
                  double dx,
                  int ny,
                  double ymin,
                  double dy,
                  int nz,
                  double zmin,
                  double dz,
                  BFMapType::enum_type atype,
                  double scale,
                  BFInterpolationStyle style,
                  bool warnIfOutside = false)
            : BFMap(filename,
                    xmin,
                    xmin + (nx - 1) * dx,
                    ymin,
                    ymin + (ny - 1) * dy,
                    zmin,
                    zmin + (nz - 1) * dz,
                    atype,
                    scale,
                    warnIfOutside),
              _nx(nx),
              _ny(ny),
              _nz(nz),
              _dx(dx),
              _dy(dy),
              _dz(dz),
              _field(_nx, _ny, _nz),
              _isDefined(_nx, _ny, _nz, false),
              _allDefined(false),
              _interpStyle(style){};

        bool getBFieldWithStatus(const CLHEP::Hep3Vector&, CLHEP::Hep3Vector&) const override;

        // Validity checker
        bool isValid(const CLHEP::Hep3Vector& point) const override;
        bool isValid(const GridPoint& ipoint) const {
            return _field.isValid(ipoint.ix, ipoint.iy, ipoint.iz);
        }

        unsigned int nx() const { return _nx; }
        unsigned int ny() const { return _ny; }
        unsigned int nz() const { return _nz; }

        double dx() const { return _dx; };
        double dy() const { return _dy; };
        double dz() const { return _dz; };

        CLHEP::Hep3Vector grid2point(unsigned ix, unsigned iy, unsigned iz) const {
            return CLHEP::Hep3Vector(_xmin + ix * _dx, _ymin + iy * _dy, _zmin + iz * _dz);
        }

        GridPoint point2grid(const CLHEP::Hep3Vector& pos) const;

        // returns vector from ipos to pos normalized to grid spacing
        CLHEP::Hep3Vector cellFraction(const CLHEP::Hep3Vector& pos, const GridPoint& ipos) const;

        // public function for getNeighbor
        bool getNeighborPointBF(const CLHEP::Hep3Vector&,
                                CLHEP::Hep3Vector neighborPoints[3],
                                CLHEP::Hep3Vector neighborBF[3][3][3]) const;

        void print(std::ostream& os) const override;

       private:
        // Grid dimensions
        unsigned int _nx, _ny, _nz;

        // Distance between points.
        double _dx, _dy, _dz;

        // Vector arrays for gridpoints and field values
        mu2e::Container3D<CLHEP::Hep3Vector> _field;
        mu2e::Container3D<bool> _isDefined;

        // If all grid points are valid then _isDefined is not needed.
        bool _allDefined;

        // Flag to flip Y component for maps that assume XZ-plane symmetry.
        bool _flipy = true;

        // method for interpolation between field grid points
        BFInterpolationStyle _interpStyle;

        // Functions used internally and by the code that populates the maps.

        // method to store the neighbors
        bool getNeighbors(int ix, int iy, int iz, CLHEP::Hep3Vector neighborsBF[3][3][3]) const;

        // Interpolator
        CLHEP::Hep3Vector interpolate(CLHEP::Hep3Vector const vec[3][3][3],
                                      const CLHEP::Hep3Vector& frac) const;

        // Polynomial fit function used by interpolator
        double gmcpoly2(double const f1d[3], double const& x) const;

        // Compute grid indices for a given point.
        std::size_t iX(double x) const { return static_cast<std::size_t>(lround((x - _xmin) / _dx)); }

        std::size_t iY(double y) const { return static_cast<std::size_t>(lround((y - _ymin) / _dy)); }

        std::size_t iZ(double z) const { return static_cast<std::size_t>(lround((z - _zmin) / _dz)); }

        bool interpolateTriLinear(const CLHEP::Hep3Vector&, CLHEP::Hep3Vector&) const;

    };

    inline BFGridMap::GridPoint BFGridMap::point2grid(const CLHEP::Hep3Vector& pos) const {
        return GridPoint(iX(pos.x()), iY(pos.y()), iZ(pos.z()));
    }

    std::ostream& operator<<(std::ostream& os, const BFGridMap::GridPoint& ipoint);

}  // end namespace mu2e

#endif /* BFieldGeom_BFGridMap_hh */
