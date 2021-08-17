#ifndef RecoDataProducts_GenVector
#define RecoDataProducts_GenVector
// typedef for cartesian vector used in reconstruction
/// root
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include <string>
using ROOT::Math::XYZVectorF;
using ROOT::Math::XYZTVectorF;
using ROOT::Math::XYZVectorD;
using XYZTVectorD = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D <double> >;
namespace mu2e {
  namespace GenVector {
    // provide a generic translation from CLHEP
    // the following is deprecated, XYZVec templated constructor now does this
    CLHEP::Hep3Vector Hep3Vec(XYZVectorF const& rvec);
    CLHEP::Hep3Vector Hep3Vec(XYZVectorD const& rvec);
    CLHEP::HepLorentzVector HepLorentzVec(XYZTVectorF const& rvec);
    CLHEP::HepLorentzVector HepLorentzVec(XYZTVectorD const& rvec);
    // z direction definition; this is missing from GenVector
    XYZVectorF const& ZDir();
  }
}
#endif
