#ifndef BFieldGeom_BFMap_hh
#define BFieldGeom_BFMap_hh
//
// Class to hold one magnetic field map. The map is defined on a regular cartesian grid.
// All field maps are given in the standard Mu2e coordinate system.
// Units are: space point in mm, field values in tesla.
//
// $Id: BFMap.hh,v 1.11 2011/05/17 15:41:35 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:35 $
//
// Original Rob Kutschke, based on work by Julie Managan and Bob Bernstein.
// Rewritten in part by Krzysztof Genser to save execution time
//

#include <iosfwd>
#include <string>
#include <vector>
#include "BFieldGeom/inc/BFMapBase.hh"
#include "BFieldGeom/inc/BFMapType.hh"
#include "BFieldGeom/inc/Container3D.hh"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {
  class BFMap : public BFMapBase {

  public:

    friend class BFieldManagerMaker;

    BFMap():
      _key(),
      _warnIfOutside(false),
      _nx(),
      _ny(),
      _nz(),
      _grid(),
      _field(),
      _isDefined(),
      _allDefined(false),
      _type(),
      _scaleFactor(1.){
    }

    BFMap( const std::string& key, BFMapType::enum_type atype, double scale, bool warnIfOutside=false):
      _key(key),
      _warnIfOutside(warnIfOutside),
      _nx(),
      _ny(),
      _nz(),
      _grid(),
      _field(),
      _isDefined(),
      _allDefined(false),
      _type(atype),
      _scaleFactor(scale){
    }

    BFMap(std::string filename, 
          int const nx, 
          int const ny, 
          int const nz,
          BFMapType::enum_type atype,
          double scale,
          bool warnIfOutside=false):
      _key(filename),
      _warnIfOutside(warnIfOutside),
      _nx(nx),
      _ny(ny),
      _nz(nz),
      _grid(_nx,_ny,_nz),
      _field(_nx,_ny,_nz),
      _isDefined(_nx,_ny,_nz,false),
      _allDefined(false),
      _type(atype),
      _scaleFactor(scale){
    };
    
    virtual ~BFMap();

    // Accessors
    bool getBFieldWithStatus(const CLHEP::Hep3Vector &, CLHEP::Hep3Vector &) const;

    // Validity checker
    bool isValid(CLHEP::Hep3Vector const& point) const;

    // Some extra checks for GMC format maps.
    bool isGMCValid(CLHEP::Hep3Vector const& point) const;

    int nx() const { return _nx; }
    int ny() const { return _ny; }
    int nz() const { return _nz; }

    double xmin() const {return _xmin;}; double xmax() const {return _xmax;};
    double ymin() const {return _ymin;}; double ymax() const {return _ymax;};
    double zmin() const {return _zmin;}; double zmax() const {return _zmax;};

    double dx() const {return _dx;}; 
    double dy() const {return _dy;};
    double dz() const {return _dz;};

    BFMapType type() const { return _type; }

    const std::string& getKey() const { return _key; };

    void print( std::ostream& os) const;

  private:

    // Filename, database key or other id information that describes
    // where this map came from.
    std::string _key;

    // If true, then print a warning message when a point is outside the region 
    // in which the map is defined; else return a field with a value of (0.,0.,0.);
    // This does happen under normal operation of G4 so we should not warn by default.
    bool _warnIfOutside;

    // Grid dimensions
    const unsigned int _nx, _ny, _nz;

    // Min and Max values.
    double _xmin, _xmax, _ymin, _ymax, _zmin, _zmax;

    // Distance between points.
    double _dx, _dy, _dz;

    // Vector arrays for gridpoints and field values
    mu2e::Container3D<CLHEP::Hep3Vector> _grid;
    mu2e::Container3D<CLHEP::Hep3Vector> _field;
    mu2e::Container3D<bool> _isDefined;

    // If all grid points are valid then _isDefined is not needed.
    bool  _allDefined;

    // GMC, G4BL or possible future types.
    BFMapType _type;

    // A scale factor applied overall.
    double _scaleFactor;

    // Functions used internally and by the code that populates the maps.

    // method to store the neighbors
    bool getNeighbors(int ix, int iy, int iz, CLHEP::Hep3Vector neighborsBF[3][3][3]) const;

    // Interpolator
    CLHEP::Hep3Vector interpolate(CLHEP::Hep3Vector const vec[3][3][3],
                                  double const frac[3]) const;

    // Polynomial fit function used by interpolator
    double gmcpoly2(double const f1d[3], double const& x) const;

    // Define the limits and step sizes for the maps.
    void setLimits ( double xmin, double xmax, double ymin, double ymax,
                     double zmin, double zmax);

    // Compute grid indices for a given point.
    std::size_t iX( double x){
      return static_cast<int>((x - _xmin)/_dx + 0.5);
    }

    std::size_t iY( double y){
      return static_cast<int>((y - _ymin)/_dy + 0.5);
    }

    std::size_t iZ( double z){
      return static_cast<int>((z - _zmin)/_dz + 0.5);
    }

  };

} // end namespace mu2e

#endif /* BFieldGeom_BFMap_hh */
