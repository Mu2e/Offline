//
//  Description  a geometric ray as a point and direction
//  original author: David Brown (LBN) 2023
//
#ifndef RecoGeom_Ray_hh
#define RecoGeom_Ray_hh
#include "Math/Vector3D.h"
#include <ostream>
using ROOT::Math::XYZVectorD;
namespace mu2e {
  namespace RecoGeom {
    struct Ray {
      // construct, making sure direciton is unit
      Ray(XYZVectorD const& dir, XYZVectorD const& start) : dir_(dir.Unit()), start_(start){}
      XYZVectorD dir_; // direction
      XYZVectorD start_; // starting position
      XYZVectorD position(double distance) const { return start_ + distance*dir_; }
    };
  }
}
std::ostream& operator <<(std::ostream& ost, mu2e::RecoGeom::Ray const& ray);
#endif
