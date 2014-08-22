#ifndef threevec_HH
#define threevec_HH
#include "Rtypes.h"
namespace mu2e
{
  struct threevec {
    Float_t _x,_y,_z;
    threevec(): _x(0.0),_y(0.0),_z(0.0) {}
    threevec(const CLHEP::Hep3Vector& vec) : _x(vec.x()),_y(vec.y()),_z(vec.z()) {}
    threevec& operator = (const CLHEP::Hep3Vector& vec){ _x =vec.x(); _y =vec.y(); _z= vec.z(); return *this; }
    void reset() { _x = _y = _z = 0.0; }
  };
}
#endif


