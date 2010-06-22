//
// Class to hold one magnetic field map. The map
// is defined on a regular cartesian grid.
//
// $Id: BFMap.cc,v 1.1 2010/06/22 16:44:25 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/06/22 16:44:25 $
//
// Original Rob Kutschke, based on work by Julie Managan and Bob Bernstein.
//

// C++ includes
#include <iostream>
#include <string>
//#include <fcntl.h>
#include <vector>
//#include <set>

// Framework includes
//#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "BFieldGeom/inc/BFMap.hh"
#include "BFieldGeom/inc/DiskRecord.hh"
#include "BFieldGeom/inc/Container3D.hh"

// Other includes
#include "CLHEP/Vector/ThreeVector.h"

using CLHEP::Hep3Vector;
using mu2e::Container3D;
using namespace std;

namespace mu2e {

  // function to determine if the point is in the map
  bool BFMap::isValid(CLHEP::Hep3Vector const& point) const{
    if (point.x() < _xmin || point.x() > _xmax) {
      cout << "Bad X: " << point << endl;
      return false;
    }
    if (point.y() < _ymin || point.y() > _ymax) {
      cout << "Bad Y: " << point << endl;
      return false;
    }
    if (point.z() < _zmin || point.z() > _zmax) {
      cout << "Bad Z: " << point << endl;
      return false;
    }
    return true;
  }

  // function to save an array of B-field values at the nearest 
  // neighbors to the point
  void BFMap::getNeighbors(Container3D<CLHEP::Hep3Vector>& neighborsBF) const {
    for (int i = 0; i != 3; ++i){
      unsigned int xindex = _ix + i - 1;
      for (int j = 0; j != 3; ++j){
        unsigned int yindex = _iy + j - 1;
        for (int k = 0; k != 3; ++k){
          unsigned int zindex = _iz + k - 1;
          neighborsBF.set(i, j, k, _field(xindex, yindex, zindex));
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
  Hep3Vector BFMap::interpolate(Container3D<CLHEP::Hep3Vector> const vec,
                                  CLHEP::Hep3Vector const frac) const {
    // Create vecs and vectors
    vector<double> x1d, y1d, z1d;
    vector<CLHEP::Hep3Vector> vecx;
    vector<CLHEP::Hep3Vector> vecxy;
    double xin = frac.x();
    double yin = frac.y();
    double zin = frac.z();

    /*  cout << "--------Starting interpolation-------- \n" << "Input: (" 
        << xin << "," << yin << "," << zin << ")" << endl;*/

    // First loop - interpolate xin
    for (int j = 0; j != 3; ++j){
      for (int k = 0; k != 3; ++k){
        for (int i = 0; i != 3; ++i){
          /*        cout << "1st loop access: (" << vec(i,j,k).x() << ","
                    << vec(i,j,k).y() << "," << vec(i,j,k).z() << ")" << endl;*/
          x1d.push_back(vec(i,j,k).x());
          y1d.push_back(vec(i,j,k).y());
          z1d.push_back(vec(i,j,k).z());
        }
        double xval = gmcpoly2(x1d,xin);
        double yval = gmcpoly2(y1d,xin);
        double zval = gmcpoly2(z1d,xin);
        /*      cout << "Xin pass: (" << xval << "," << yval << "," << zval
                << ")" << endl;*/
        vecx.push_back(CLHEP::Hep3Vector(xval,yval,zval));
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
      double xval = gmcpoly2(x1d,yin);
      double yval = gmcpoly2(y1d,yin);
      double zval = gmcpoly2(z1d,yin);
      vecxy.push_back(CLHEP::Hep3Vector(xval,yval,zval));
    }

    // Third loop - interpolate zin
    for (int k = 0; k != 3; ++k){
      x1d[k] = vecxy[k].x();
      y1d[k] = vecxy[k].y();
      z1d[k] = vecxy[k].z();
    }
    double xval = gmcpoly2(x1d,zin);
    double yval = gmcpoly2(y1d,zin);
    double zval = gmcpoly2(z1d,zin);

    // Return BField 3Vector
    return CLHEP::Hep3Vector(xval,yval,zval);
  }

  // Standard Lagrange formula for 2nd order polynomial fit of 
  // univariate function
  double BFMap::gmcpoly2(vector<double> const& f1d, double const& x) const {

    double x0(0.),x1(1.),x2(2.);
    double y0 = f1d[0];
    double y1 = f1d[1];
    double y2 = f1d[2];

    double fout = y0*(x-x1)*(x-x2)/((x0-x1)*(x0-x2)) +
      y1*(x-x0)*(x-x2)/((x1-x0)*(x1-x2)) +
      y2*(x-x0)*(x-x1)/((x2-x0)*(x2-x1));
    return fout;
  } 

  // Function to return the BField for any point
  Hep3Vector BFMap::getBField(CLHEP::Hep3Vector const& testpoint ) const{

    // Allow y-symmetry (grid is only defined for y > 0);
    int sign(1);
    CLHEP::Hep3Vector point(testpoint.x(), testpoint.y(), testpoint.z());
    if (testpoint.y() < _origin.y()){
      sign = -1;
      double y = -(testpoint.y() - _origin.y()) + _origin.y();
      point.setY(y);
    }
    //    cout << "Point: " << point << endl;
    //cout << "Point: " << point+CLHEP::Hep3Vector(0.,-7350.,-4000.) << endl;
    //cout << "Origin: " << _origin << endl;

    // Check validity.  Throw or return a zero field.
    if ( !isValid(point) ){
      if ( _throwIfOutside ){
        throw cms::Exception("GEOM")
          << "Point is outside of the valid region of the map: " << _key << "\n"
          << "Point in input coordinates: " << testpoint << "\n";
      }
      return CLHEP::Hep3Vector(0.,0.,0.);
    }

    // Get the indices of the nearest grid point
    _ix = static_cast<int>((point.x() - _xmin)/_dx + 0.5);
    _iy = static_cast<int>((point.y() - _ymin)/_dy + 0.5);
    _iz = static_cast<int>((point.z() - _zmin)/_dz + 0.5);
    /*
    cout << "Initial ix,y,z: " 
         << _ix << " "
         << _iy << " "
         << _iz << " "
         << endl;
    */

    // Correct for edge points by moving their NGPt just inside the edge
    if (_ix ==     0){++_ix;}
    if (_ix == _nx-1){--_ix;}
    if (_iy ==     0){++_iy;}
    if (_iy == _ny-1){--_iy;}
    if (_iz ==     0){++_iz;}
    if (_iz == _nz-1){--_iz;}
    /*
    cout << "Final ix,y,z:  " 
         << _ix << " "
         << _iy << " "
         << _iz << " "
         << endl;
    */
  
    /*  cout << "Nearest Point: (" << _grid(_ix,_iy,_iz).x() << ","
        << _grid(_ix,_iy,_iz).y() << "," << _grid(_ix,_iy,_iz).z() 
        << ")\nIndices set to: " << _ix << " " << _iy << " " << _iz 
        << "\nField: (" << _field(_ix,_iy,_iz).x() << ","
        << _field(_ix,_iy,_iz).y() << "," << _field(_ix,_iy,_iz).z() << ")"
        << endl;
    */ 
    // Get the BField values of the nearest grid neighbors to the point
    Container3D<CLHEP::Hep3Vector> neighborsBF(3,3,3);
    getNeighbors(neighborsBF);

    // Set up fractional grid point. The interpolater treats the neighbors array
    // as spanning from (0,0,0) to (2,2,2), so we find the location of the point
    // within this frame.

    unsigned int xindex = _ix-1; 
    unsigned int yindex = _iy-1; 
    unsigned int zindex = _iz-1;
    double xFrac = (point.x() - _grid(xindex,yindex,zindex).x())/_dx;
    double yFrac = (point.y() - _grid(xindex,yindex,zindex).y())/_dy;
    double zFrac = (point.z() - _grid(xindex,yindex,zindex).z())/_dz;
    /*
    cout << "Frac: "
         << xFrac << " "
         << yFrac << " "
         << zFrac << " "
         << endl;
    */

    CLHEP::Hep3Vector frac = CLHEP::Hep3Vector(xFrac,yFrac,zFrac);
    /*  cout << "Frac: (" << frac.x() << "," << frac.y() << ","
        << frac.z() << ")" << endl;*/

    // Run the interpolater
    CLHEP::Hep3Vector testBF = interpolate(neighborsBF,frac);
    //    cout << "testBF: " << testBF << endl;

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
    cout << "Size:       " << _dx << " " << _dy << " " << _dz << endl;

    if ( _throwIfOutside ){
      cout << "Will throw if outside of the valid region." << endl;
    } else{
      cout << "Will return (0.,0.,0.) if outside of the valid region." << endl;
    }
  }

} // end namespace mu2e
