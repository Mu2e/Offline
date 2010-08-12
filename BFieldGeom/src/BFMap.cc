//
// Class to hold one magnetic field map. The map
// is defined on a regular cartesian grid.
//
// $Id: BFMap.cc,v 1.4 2010/08/12 19:39:58 genser Exp $
// $Author: genser $
// $Date: 2010/08/12 19:39:58 $
//
// Original Rob Kutschke, based on work by Julie Managan and Bob Bernstein.
// Rewritten in part by Krzysztof Genser to correct mistake pointed by RB and to save execution time
//

// C++ includes
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Mu2e includes
#include "BFieldGeom/inc/BFMap.hh"
#include "BFieldGeom/inc/DiskRecord.hh"
#include "BFieldGeom/inc/Container3D.hh"

// Other includes
#include "CLHEP/Vector/ThreeVector.h"

//using CLHEP::Hep3Vector;
//using mu2e::Container3D;
using namespace std;

namespace mu2e {

  // function to determine if the point is in the map
  bool BFMap::isValid(CLHEP::Hep3Vector const& point) const{
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

  // Populate a 3x3x3 array with the field values centered on the grid point
  // that is the closest to the requested point.
  void BFMap::getNeighbors( int ix, int iy, int iz, 
                            CLHEP::Hep3Vector neighborsBF[3][3][3]) const {
    for (int i = 0; i != 3; ++i){
      unsigned int xindex = ix + i - 1;
      for (int j = 0; j != 3; ++j){
        unsigned int yindex = iy + j - 1;
        for (int k = 0; k != 3; ++k){
          unsigned int zindex = iz + k - 1;
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
    return;
  }

  // Function to interpolate the BField value at the point from the values
  // at in the neighbor grid. 
  CLHEP::Hep3Vector BFMap::interpolate(CLHEP::Hep3Vector const vec[3][3][3],
                                       double const frac[3]) const {

    // Create vecs and vectors
    double x1d[3], y1d[3], z1d[3];
    static CLHEP::Hep3Vector vecx[9];
    static CLHEP::Hep3Vector vecxy[3];
    double xin = frac[0];
    double yin = frac[1];
    double zin = frac[2];

    /*  cout << "--------Starting interpolation-------- \n" << "Input: (" 
        << xin << "," << yin << "," << zin << ")" << endl;*/

    // First loop - interpolate xin
    for (int j = 0; j != 3; ++j){
      for (int k = 0; k != 3; ++k){
        for (int i = 0; i != 3; ++i){
          /*        cout << "1st loop access: (" << vec[i][j][k].x() << ","
                    << vec[i][j][k].y() << "," << vec[i][j][k].z() << ")" << endl;*/
          x1d[i] = (vec[i][j][k]).x();
          y1d[i] = (vec[i][j][k]).y();
          z1d[i] = (vec[i][j][k]).z();
        }
        vecx[j*3 + k] = CLHEP::Hep3Vector(gmcpoly2(x1d,xin),
                                          gmcpoly2(y1d,xin),
                                          gmcpoly2(z1d,xin));
      }
    }

    // Second loop - interpolate yin
    for (int k = 0; k != 3; ++k){
      for (int j = 0; j != 3; ++j){
        int index = j*3 + k;
        x1d[j] = vecx[index].x();
        y1d[j] = vecx[index].y();
        z1d[j] = vecx[index].z();
      }
      vecxy[k] = CLHEP::Hep3Vector(gmcpoly2(x1d,yin),
                                   gmcpoly2(y1d,yin),
                                   gmcpoly2(z1d,yin));
    }

    // Third loop - interpolate zin
    for (int k = 0; k != 3; ++k){
      x1d[k] = vecxy[k].x();
      y1d[k] = vecxy[k].y();
      z1d[k] = vecxy[k].z();
    }

    // Return BField 3Vector
    return CLHEP::Hep3Vector(gmcpoly2(x1d,zin),
                             gmcpoly2(y1d,zin),
                             gmcpoly2(z1d,zin));
  }

  // Standard Lagrange formula for 2nd order polynomial fit of 
  // univariate function
  double BFMap::gmcpoly2(double const f1d[3], double const& x) const {

    static const double  x0(0.),x1(1.),x2(2.);

    return f1d[0]*(x-x1)*(x-x2)/((x0-x1)*(x0-x2)) +
      f1d[1]*(x-x0)*(x-x2)/((x1-x0)*(x1-x2)) +
      f1d[2]*(x-x0)*(x-x1)/((x2-x0)*(x2-x1));

  } 

  // Function to return the BField for any point
  CLHEP::Hep3Vector BFMap::getBField(CLHEP::Hep3Vector const& testpoint ) const{

    static const bool dflag = false;

    // Allow y-symmetry (grid is only defined for y > 0);
    int sign(1);
    CLHEP::Hep3Vector point(testpoint.x(), testpoint.y(), testpoint.z());
    if (testpoint.y() < _origin.y()){
      sign = -1;
      double y = -(testpoint.y() - _origin.y()) + _origin.y();
      point.setY(y);
    }
    
    dflag && cout 
      << "Testpoint:          " << testpoint << endl
      << "Point:              " << point << " "
      << "Origin:             " << _origin << endl;


//     LogDebug("BFIELD")
//       << "Testpoint:          " << testpoint << " "
//       << "Point:              " << point << " "
//       << "Origin:             " << _origin;

    // Check validity.  Return a zero field and optionally print a warning.
    if ( !isValid(point) ){
      if ( _warnIfOutside ){
        edm::LogWarning("GEOM")
          << "Point is outside of the valid region of the map: " << _key << "\n"
          << "Point in input coordinates: " << testpoint << "\n";
      }
      return CLHEP::Hep3Vector(0.,0.,0.);
    }

    // Get the indices of the nearest grid point
    int ix = static_cast<int>((point.x() - _xmin)/_dx + 0.5);
    int iy = static_cast<int>((point.y() - _ymin)/_dy + 0.5);
    int iz = static_cast<int>((point.z() - _zmin)/_dz + 0.5);
    
    dflag && cout << "Initial ix,iy,iz " << endl
         << "& dx, dy, dy:    " 
         << setw(4) << ix << " "  
         << setw(4) << iy << " " 
         << setw(4) << iz << " " << endl
         << setw(4) << _dx << " "
         << setw(4) << _dy << " "
         << setw(4) << _dz << " "
         << endl;
    
    // Correct for edge points by moving their NGPt just inside the edge
    if (ix ==     0){++ix;}
    if (ix == _nx-1){--ix;}
    if (iy ==     0){++iy;}
    if (iy == _ny-1){--iy;}
    if (iz ==     0){++iz;}
    if (iz == _nz-1){--iz;}
    
    dflag && cout << "Final ix,iy,iz:  " 
         << setw(4) << ix << " "
         << setw(4) << iy << " "
         << setw(4) << iz << " "
         << endl;
    
  
    dflag && cout << "Nearest Point:   " << _grid(ix,iy,iz) << endl
         << "Indices set to:  " << setw(4) << ix << " " << setw(4) << iy << " " << setw(4) << iz << endl
         << "Field:              " << _field(ix,iy,iz) << endl;

    // Get the BField values of the nearest grid neighbors to the point
    static CLHEP::Hep3Vector neighborsBF[3][3][3];
    getNeighbors(ix, iy, iz, neighborsBF);

    // Set up fractional grid point. The interpolator treats the neighbors array
    // as spanning from (0,0,0) to (2,2,2), so we find the location of the point
    // within this frame.

    // klg it is done to be consistent with the getNeighbors point selection 
    // klg doing it on the unit grid saves passing the _grid neighbor points to the interpolator

    int xindex = ix-1;
    int yindex = iy-1; 
    int zindex = iz-1;

    dflag && cout << "Used   x, y, z:  " 
         << setw(4) << xindex << " "
         << setw(4) << yindex << " "
         << setw(4) << zindex << " "
         << endl;

    dflag && cout << "Used Point:      " << _grid(xindex,yindex,zindex) << endl;

    double frac[] = {
      (point.x() - _grid(xindex,yindex,zindex).x())/_dx,
      (point.y() - _grid(xindex,yindex,zindex).y())/_dy,
      (point.z() - _grid(xindex,yindex,zindex).z())/_dz 
    };

    // Run the interpolator
    CLHEP::Hep3Vector testBF = interpolate(neighborsBF,frac);
    dflag && cout << "Interpolated Field: " << testBF << endl;

    // Reassign y sign
    if (sign == -1){
      testBF.setY(-testBF.y());
    }
    return testBF;
  }

  // Called by BFieldManagerMaker.
  void BFMap::setLimits ( double xmin, double xmax, double ymin, double ymax,
                          double zmin, double zmax){
    _xmin = xmin;
    _xmax = xmax;
    _ymin = ymin;
    _ymax = ymax;
    _zmin = zmin;
    _zmax = zmax;

    // Set the distance between gridpoints
    _dx = (_xmax - _xmin)/(_nx-1.);
    _dy = (_ymax - _ymin)/(_ny-1.);
    _dz = (_zmax - _zmin)/(_nz-1.);

  }

  void BFMap::print( std::ostream& os) const{
    os << "Magnetic Field Info: " << _key << endl;

    cout << "Dimensions: " << _nx << " " << _ny << " " << _nz << endl;
    cout << "Range X:    " << _xmin << " : " << _xmax << "  Middle: " << (_xmin+_xmax)/2. << endl;
    cout << "Range Y:    " << _ymin << " : " << _ymax << "  Middle: " << (_ymin+_ymax)/2. << endl;
    cout << "Range Z:    " << _zmin << " : " << _zmax << "  Middle: " << (_zmin+_zmax)/2. << endl;
    cout << "Distance:       " << _dx << " " << _dy << " " << _dz << endl;

    cout << "Field at the edges: " 
         << _field(    0,     0,     0) << ", "
         << _field(_nx-1,     0,     0) << ", "
         << _field(    0, _ny-1,     0) << ", "
         << _field(    0,     0, _nz-1) << ", "
         << _field(_nx-1, _ny-1,     0) << ", "
         << _field(_nx-1, _ny-1, _nz-1) << endl;

    cout << "Field in the middle: "
         << _field(_nx/2, _ny/2, _nz/2) << endl;

    if ( _warnIfOutside ){
      cout << "Will warn if outside of the valid region." << endl;
    } else{
      cout << "Will not warn if outside of the valid region." << endl;
    }
  }

} // end namespace mu2e
