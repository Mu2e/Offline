#ifndef BFWCONT_HH
#define BFWCONT_HH

#include <string>
#include <vector>
#include "Mu2eG4/inc/Record.hh"
#include "Mu2eG4/inc/Container3D.hh"
#include "CLHEP/Vector/ThreeVector.h"

class BFwcont { //public BField

public:
  BFwcont(std::string filename, CLHEP::Hep3Vector const& origin,
	  int const nx, int const ny, int const nz ): 
    _nx(nx),
    _ny(ny),
    _nz(nz),
    _origin(origin),
    _grid(_nx,_ny,_nz),
    _field(_nx,_ny,_nz){

    readmap(filename);
  };

  ~BFwcont(){
  };

  // Accessors
  CLHEP::Hep3Vector GetBField(CLHEP::Hep3Vector const& point);

  // function to read in the map
  void readmap(std::string filename);

  double xmin() const {return _xmin;}; double xmax() const {return _xmax;};
  double ymin() const {return _ymin;}; double ymax() const {return _ymax;};
  double zmin() const {return _zmin;}; double zmax() const {return _zmax;};

  double dx() const {return _dx;}; 
  double dy() const {return _dy;};
  double dz() const {return _dz;};

private:
  // Grid dimensions
  const unsigned int _nx, _ny, _nz;

  // Min and Max values (set in readmap)
  double _xmin, _xmax, _ymin, _ymax, _zmin, _zmax;

  // Distance between points (set in readmap)
  double _dx, _dy, _dz;

  // Validity switch and indices of nearest gridpoint
  bool _inMap;
  unsigned int _ix, _iy, _iz;

  // Origin from external setup
  CLHEP::Hep3Vector _origin;

  // Vector arrays for gridpoints and field values
  mu2e::Container3D<CLHEP::Hep3Vector> _grid; 
  mu2e::Container3D<CLHEP::Hep3Vector> _field;

  // Validity checker
  bool isValid(CLHEP::Hep3Vector const& point);

  // method to store the neighbors
  void getNeighbors(mu2e::Container3D<CLHEP::Hep3Vector>& neighborsBF) const;

  // Interpolater
  CLHEP::Hep3Vector interpolate(mu2e::Container3D<CLHEP::Hep3Vector> const vec,
				CLHEP::Hep3Vector const frac) const;

  // Polynomial fit function used by interpolater
  double gmcpoly2(std::vector<double> const& f1d, double const& x) const;

};

#endif
