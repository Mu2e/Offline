#ifndef RecoDataProducts_XYZVec_hh
#define RecoDataProducts_XYZVec_hh
// typedef for cartesian vector used in reconstruction
/// root
#include "Math/Vector3D.h"
#include "CLHEP/Vector/ThreeVector.h"
typedef ROOT::Math::XYZVectorF  XYZVec;
// provide a generic translation from CLHEP
XYZVec toXYZVec(CLHEP::Hep3Vector const& cvec) { return
  XYZVec((float)cvec.x(),(float)cvec.y(),(float)cvec.z()); }
CLHEP::Hep3Vector toCVec(XYZVec const& rvec) { return
  CLHEP::Hep3Vector((double)rvec.x(),(double)rvec.y(),(double)rvec.z()); }
// z direction definition; this is missing from GenVector
namespace Geom {
  XYZVec const& ZDir() {
    static XYZVec _zdir(0.0,0.0,1.0);
    return _zdir;
  }
}
#endif
