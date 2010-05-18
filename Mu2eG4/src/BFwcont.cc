//
// Class to access a MECO style magnetic field map.
//
// $Id: BFwcont.cc,v 1.4 2010/05/18 22:33:49 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/05/18 22:33:49 $
//
// Original author Julie Managan and Bob Bernstein

#include <iostream>
#include <string>
#include <fcntl.h>
#include <vector>
#include <set>
#include "Mu2eG4/inc/BFwcont.hh"
#include "Mu2eG4/inc/DiskRecord.hh"
#include "Mu2eG4/inc/Record.hh"
#include "Mu2eG4/inc/Container3D.hh"
#include "CLHEP/Vector/ThreeVector.h"

using CLHEP::Hep3Vector;
using mu2e::Container3D;
using namespace std;

namespace mu2e {

  // Function to read in the field map and set variables
  void BFwcont::readmap(string fieldmap){

    // Open the input file.
    int fd = open( fieldmap.c_str(), O_RDONLY );
    cout << "Open status: " << fd << endl;
    if ( !fd ) exit(-1);
  
    // data holders
    vector<Record> data;
    set<double> X, Y, Z;

    // Working buffer.
    DiskRecord r;

    // Read in the file line by line into the Diskrecord
    int nrec(0);
    while (true){
      ssize_t s = read( fd, &r, sizeof(DiskRecord) );
      if ( !s ) {
        cout << "EOF after reading record: " << nrec << endl;
        break;
      }
      ++nrec;

      // unit conversion
      r.x *= 10; r.y *= 10; r.z *= 10;
      r.bx /= 10; r.by /= 10; r.bz /= 10;

      // re-centering on (0,0,11000)
      r.x += _origin.x(); r.y += _origin.y(); r.z += _origin.z();

      // The one check I can do.
      if ( r.head != r.tail ){
        cout << "Mismatched head and tail at record: " << nrec << endl;
        exit(-1);
      }

      // Save the record.
      data.push_back( Record(r) );
      X.insert(r.x);
      Y.insert(r.y);
      Z.insert(r.z);
    }

    vector<double> vX(X.begin(),X.end());
    vector<double> vY(Y.begin(),Y.end());
    vector<double> vZ(Z.begin(),Z.end());

    // Set the min and max values
    _xmin = vX.front();  _xmax = vX.back();
    _ymin = vY.front();  _ymax = vY.back();
    _zmin = vZ.front();  _zmax = vZ.back();
    cout << "X: " << _xmin << " to " << _xmax << "\nY: " << _ymin 
         << " to " << _ymax << "\nZ: " << _zmin << " to " << _zmax 
         << endl;

    // Set the distance between gridpoints
    _dx = (_xmax - _xmin)/(_nx-1.);
    _dy = (_ymax - _ymin)/(_ny-1.);
    _dz = (_zmax - _zmin)/(_nz-1.);

    // Store grid and field 3D arrays
    for (vector<Record>::const_iterator i=data.begin(), e=data.end();
         i != e; ++i){
      Record const& r(*i);
      unsigned int xIndex = static_cast<int>((r.x - _xmin)/_dx + 0.5);
      unsigned int yIndex = static_cast<int>((r.y - _ymin)/_dy + 0.5);
      unsigned int zIndex = static_cast<int>((r.z - _zmin)/_dz + 0.5);

      _grid.set(xIndex, yIndex, zIndex, CLHEP::Hep3Vector(r.x,r.y,r.z));
      _field.set(xIndex, yIndex, zIndex, CLHEP::Hep3Vector(r.bx,r.by,r.bz));

    }

    return;
  }

  // function to determine if the point is in the map
  bool BFwcont::isValid(CLHEP::Hep3Vector const& point){
    if (point.x() < _xmin || point.x() > _xmax) {
      cout << "Bad X: " << point.x() << endl;
      return false;
    }
    if (point.y() < _ymin || point.y() > _ymax) {
      cout << "Bad Y: " << point.y() << endl;
      return false;
    }
    if (point.z() < _zmin || point.z() > _zmax) {
      cout << "Bad Z: " << point.z() << endl;
      return false;
    }
    return true;
  }

  // function to save an array of B-field values at the nearest 
  // neighbors to the point
  void BFwcont::getNeighbors(Container3D<CLHEP::Hep3Vector>& neighborsBF) const {
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
  Hep3Vector BFwcont::interpolate(Container3D<CLHEP::Hep3Vector> const vec,
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
  double BFwcont::gmcpoly2(vector<double> const& f1d, double const& x) const {

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
  Hep3Vector BFwcont::GetBField(CLHEP::Hep3Vector const& testpoint ){

    // Allow y-symmetry (grid is only defined for y > 0);
    int sign(1);
    CLHEP::Hep3Vector point(testpoint.x(), testpoint.y(), testpoint.z());
    if (testpoint.y() < _origin.y()){
      sign = -1;
      double y = -(testpoint.y() - _origin.y()) + _origin.y();
      point.setY(y);
    }

    // check validity
    _inMap = isValid(point);
    if (! _inMap){ 
      return CLHEP::Hep3Vector(0,0,0);
    }

    // Get the indices of the nearest grid point
    _ix = static_cast<int>((point.x() - _xmin)/_dx + 0.5);
    _iy = static_cast<int>((point.y() - _ymin)/_dy + 0.5);
    _iz = static_cast<int>((point.z() - _zmin)/_dz + 0.5);

    // Correct for edge points by moving their NGPt just inside the edge
    if (_ix ==     0){++_ix;}
    if (_ix == _nx-1){--_ix;}
    if (_iy ==     0){++_iy;}
    if (_iy == _ny-1){--_iy;}
    if (_iz ==     0){++_iz;} 
    if (_iz == _nz-1){--_iz;}
  
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

    CLHEP::Hep3Vector frac = CLHEP::Hep3Vector(xFrac,yFrac,zFrac);
    /*  cout << "Frac: (" << frac.x() << "," << frac.y() << ","
        << frac.z() << ")" << endl;*/

    // Run the interpolater
    CLHEP::Hep3Vector testBF = interpolate(neighborsBF,frac);

    // Reassign y sign
    if (sign == -1){
      testBF.setY(-testBF.y());
    }

    return testBF;
  }

} // end namespace mu2e
