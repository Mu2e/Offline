//
//  Description of 2D points with covariance matrices
//  This class works in SINGLE PRECISION for speed, so has limited accuracy
//  Original author: Dave Brown (LBNL) 4/7/2023
//
#ifndef GeneralUtilities_TwoDPoint_hh
#define GeneralUtilities_TwoDPoint_hh
#include "Math/Vector3D.h"
#include "Math/Vector2D.h"
#include "Math/SVector.h"
#include "Math/SMatrix.h"
namespace mu2e {
  class TwoDPoint {
    public:
      using VEC3 = ROOT::Math::XYZVectorF; // spatial only 3D vector
      using VEC2 = ROOT::Math::XYVectorF; // spatial only 2D vector
      using SVEC = ROOT::Math::SVector<float,2>; // 2D algebraic vector
      using SMAT = ROOT::Math::SMatrix<float,2,2,ROOT::Math::MatRepSym<float,2>>; // 2D symmetric matrix
      TwoDPoint() {}
      TwoDPoint(SVEC const& pos, SMAT const& cov) : pos_(pos), cov_(cov) {} // explicit construct
      TwoDPoint(VEC2 const& pos,float ucos, float uvar, float vvar); // construct from a point, semi-major direction cosine to the X axis, and semi-major and minor variances
      TwoDPoint(VEC3 const& pos,VEC3 const& udir, float uvar, float vvar); // construct from a point, semi-major direction, and semi-major and minor variances.  This is used for ComboHit
      // accessors
      auto const& pos() const { return pos_; }
      auto const& cov() const { return cov_; }
      SMAT weight(float intrinsicvar=0.0) const; // invert the covariance to create the weight matrix.  Optionaly add a constant variance
      VEC2 pos2() const { return VEC2(pos_[0],pos_[1]); } // physical vector
      VEC3 pos3() const { return VEC3(pos_[0],pos_[1],0.0); } // physical vector; Z coordinate is meaningless
    private:
      SVEC pos_; // position in XY cartesian coordinates represented as an algebraic vector
      SMAT cov_; // covariance matrix in XY cartesian coordinates
  };
}
#endif
