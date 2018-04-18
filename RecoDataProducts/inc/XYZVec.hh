#ifndef RecoDataProducts_XYZVec_hh
#define RecoDataProducts_XYZVec_hh
// typedef for cartesian vector used in reconstruction
/// root
#include "Math/Vector3D.h"
#include "CLHEP/Vector/ThreeVector.h"
typedef ROOT::Math::XYZVectorF  XYZVec;
namespace Geom {
// provide a generic translation from CLHEP
  XYZVec toXYZVec(CLHEP::Hep3Vector const& cvec);
  CLHEP::Hep3Vector Hep3Vec(XYZVec const& rvec);
  // z direction definition; this is missing from GenVector
  XYZVec const& ZDir();
}
#endif
