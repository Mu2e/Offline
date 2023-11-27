//
//  Description of linear track segment within a stereo hit
//  Original author: Dave Brown (LBNL) 11/26/2023
//
#ifndef TrkHitReco_StereoLine_hh
#define TrkHitReco_StereoLine_hh
#include "Math/Vector3D.h"
#include "Math/SVector.h"
#include "Math/SMatrix.h"
#include <iostream>
namespace mu2e {
  struct StereoLine {
    enum PARS {px=0, py, rx, ry}; // parameter indices
    using VEC3 = ROOT::Math::XYZVectorF; // spatial 3D vector
    using SVEC = ROOT::Math::SVector<double,4>; // algebraic vector; first 2 rows are position, last 2 are direction (dX/dZ, dY/dZ)
    using SMAT = ROOT::Math::SMatrix<double,4,4,ROOT::Math::MatRepSym<double,4>>; // matrix for weights/covariance of above
    // accessors
    auto const& pars() const { return pars_; }
    auto const& cov() const { return cov_; }
    double chisq() const { return chisq_; }
    unsigned ndof() const { return ndof_; }
    double z0() const { return z0_; }
    VEC3 pos(float zval=0.0) const; // point at a given Z position
    VEC3 dir() const; // constant direction
    // data payload
    SVEC pars_; //  parameters; Px, Py, Rx, Ry
    SMAT cov_; // covariance of the above
    float z0_; // reference z position
    double chisq_; // chisquared
    unsigned ndof_; // # of degress of freedom
  };
}
std::ostream& operator << (std::ostream& ost, mu2e::StereoLine const& sline);
#endif
