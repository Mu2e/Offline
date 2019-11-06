//
// Given a vector and a volume, gives the intersections between the two
//
//
// _original author Stefano Roberto Soleti
//

#include <iostream>
#include "Mu2eUtilities/inc/VectorVolume.hh"

using namespace std;

namespace mu2e {


  VectorVolume::VectorVolume(CLHEP::Hep3Vector const& orig,
                             CLHEP::Hep3Vector const& dir,
                             float xMin, float xMax,
                             float yMin, float yMax,
                             float zMin, float zMax):
                             _orig(orig),
                             _dir(dir),
                             _xMin(xMin),
                             _xMax(xMax),
                             _yMin(yMin),
                             _yMax(yMax),
                             _zMin(zMin),
                             _zMax(zMax)
  {
  }

  VectorVolume::~VectorVolume(){
  }

  bool VectorVolume::pointInBox(float x, float y, float x0, float y0,
                                        float x1, float y1)
  {
    bool ret = false;
    if ((x >= x0) && (x <= x1) && (y >= y0) && (y <= y1))
    {
      ret = true;
    }
    return ret;
  }

  float VectorVolume::distance(const CLHEP::Hep3Vector &u, const CLHEP::Hep3Vector &v)
  {
    return safeSqrt((u.x() - v.x()) * (u.x() - v.x()) +
                    (u.y() - v.y()) * (u.y() - v.y()) +
                    (u.z() - v.z()) * (u.z() - v.z()));
  }

  void VectorVolume::calIntersections(std::vector<CLHEP::Hep3Vector> & intersections)
  {
    // roof: _targetBox_yMax, _targetBox_xMin, _targetBox_xMax, _targetBox_zMin, _targetBox_zMax
    // skip projection if the particle goes parallely to the plane
    if (_dir.y() != 0.)
    {
      const float t = (_yMax - _orig.y()) / _dir.y();
      const float x1 = _dir.x() * t + _orig.x();
      const float z1 = _dir.z() * t + _orig.z();
      // std::cout << "x1 " << x1 << " " << z1 << std::endl;

      if (pointInBox(x1, z1, _xMin, _zMin, _xMax, _zMax))
      {
        intersections.push_back(CLHEP::Hep3Vector(x1, _yMax, z1));
      }
    }

    //east: _zMin, _xMin, _xMax, _yMin, _yMax
    if (_dir.z() != 0.)
    {
      const float t = (_zMin - _orig.z()) / _dir.z();
      const float x1 = _dir.x() * t + _orig.x();
      const float y1 = _dir.y() * t + _orig.y();
      if (pointInBox(x1, y1, _xMin, _yMin, _xMax, _yMax))
      {
        intersections.push_back(CLHEP::Hep3Vector(x1, y1, _zMin));
      }
    }

    //west: _zMax, _xMin, _xMax, _yMin, _yMax
    if (_dir.z() != 0.)
    {
      const float t = (_zMax - _orig.z()) / _dir.z();
      const float x1 = _dir.x() * t + _orig.x();
      const float y1 = _dir.y() * t + _orig.y();
      if (pointInBox(x1, y1, _xMin, _yMin, _xMax, _yMax))
      {
        intersections.push_back(CLHEP::Hep3Vector(x1, y1, _zMax));
      }
    }

    //south: _xMin, _yMin, _yMax, _zMin, _zMax
    if (_dir.x() != 0.)
    {
      const float t = (_xMin - _orig.x()) / _dir.x();
      const float z1 = _dir.z() * t + _orig.z();
      const float y1 = _dir.y() * t + _orig.y();
      if (pointInBox(z1, y1, _zMin, _yMin, _zMax, _yMax))
      {
        intersections.push_back(CLHEP::Hep3Vector(_xMin, y1, z1));
      }
    }

    //north: _xMax, _yMin, _yMax, _zMin, _zMax
    if (_dir.x() != 0.)
    {
      const float t = (_xMax - _orig.x()) / _dir.x();
      const float z1 = _dir.z() * t + _orig.z();
      const float y1 = _dir.y() * t + _orig.y();
      if (pointInBox(z1, y1, _zMin, _yMin, _zMax, _yMax))
      {
        intersections.push_back(CLHEP::Hep3Vector(_xMax, y1, z1));
      }
    }
  }

} // end namespace mu2e
