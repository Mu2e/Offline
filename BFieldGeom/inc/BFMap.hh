#ifndef BFMAP_HH
#define BFMAP_HH
//
// Class to hold one magnetic field map. The map is defined on a regular cartesian grid.
// All field maps are given in the standard Mu2e coordinate system.
// Units are: space point in mm, field values in tesla.
//
// $Id: BFMap.hh,v 1.1 2010/06/22 16:44:25 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/06/22 16:44:25 $
//
// Original Rob Kutschke, based on work by Julie Managan and Bob Bernstein.
//

#include <iosfwd>
#include <string>
#include <vector>
#include "BFieldGeom/inc/Container3D.hh"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {
  class BFMap {

  public:

    friend class BFieldManagerMaker;

    BFMap():
      _key(),
      _throwIfOutside(false),
      _nx(),
      _ny(),
      _nz(),
      _origin(),
      _grid(),
      _field(){
    }

    BFMap( const std::string& key, bool throwIfOutside=false):
      _key(key),
      _throwIfOutside(throwIfOutside),
      _nx(),
      _ny(),
      _nz(),
      _origin(),
      _grid(),
      _field(){
    }

    BFMap(std::string filename, 
          CLHEP::Hep3Vector const& origin,
          int const nx, 
          int const ny, 
          int const nz,
          bool throwIfOutside=false):
      _key(filename),
      _throwIfOutside(throwIfOutside),
      _nx(nx),
      _ny(ny),
      _nz(nz),
      _origin(origin),
      _grid(_nx,_ny,_nz),
      _field(_nx,_ny,_nz){
    };
    
    ~BFMap(){};

    // Accessors
    CLHEP::Hep3Vector getBField(CLHEP::Hep3Vector const& point) const;

    int nx() const { return _nx; }
    int ny() const { return _ny; }
    int nz() const { return _nz; }

    double xmin() const {return _xmin;}; double xmax() const {return _xmax;};
    double ymin() const {return _ymin;}; double ymax() const {return _ymax;};
    double zmin() const {return _zmin;}; double zmax() const {return _zmax;};

    double dx() const {return _dx;}; 
    double dy() const {return _dy;};
    double dz() const {return _dz;};

    const std::string& getKey() const { return _key; };

    void print( std::ostream& os) const;

  private:

    // Filename, database key or other id information that describes
    // where this map came from.
    std::string _key;

    // If true, then throw when a point is outside the region in which the 
    // map is defined; else return a field with a value of (0.,0.,0.);
    bool _throwIfOutside;

    // Grid dimensions
    const unsigned int _nx, _ny, _nz;

    // Min and Max values.
    double _xmin, _xmax, _ymin, _ymax, _zmin, _zmax;

    // Distance between points.
    double _dx, _dy, _dz;

    // Indices of nearest gridpoint
    // These should be local variables of the one member function
    // and passed as arguments to the other.
    mutable unsigned int _ix, _iy, _iz;

    // Origin from external setup
    CLHEP::Hep3Vector _origin;

    // Vector arrays for gridpoints and field values
    mu2e::Container3D<CLHEP::Hep3Vector> _grid;
    mu2e::Container3D<CLHEP::Hep3Vector> _field;

    // Functions used internally and by the code that populates the maps.

    // Validity checker
    bool isValid(CLHEP::Hep3Vector const& point) const;

    // method to store the neighbors
    void getNeighbors(mu2e::Container3D<CLHEP::Hep3Vector>& neighborsBF) const;

    // Interpolater
    CLHEP::Hep3Vector interpolate(mu2e::Container3D<CLHEP::Hep3Vector> const vec,
                                  CLHEP::Hep3Vector const frac) const;

    // Polynomial fit function used by interpolater
    double gmcpoly2(std::vector<double> const& f1d, double const& x) const;

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

#endif
