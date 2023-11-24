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
#include <iostream>
namespace mu2e {
  struct UVVar {
    using VEC2 = ROOT::Math::XYVectorF; // spatial only 2D vector
    VEC2 const& udir() const { return udir_;}
    VEC2 vdir() const { return VEC2(-udir_.Y(),udir_.X()); }
    float const& uvar() const { return uvar_; }
    float const& vvar() const { return vvar_; }
    float ures() const { return sqrt( uvar_); }
    float vres() const { return sqrt( vvar_); }
    VEC2 udir_; // u direction
    float uvar_, vvar_; // variance along u, v
  };
  class TwoDPoint {
    public:
      using VEC3 = ROOT::Math::XYZVectorF; // spatial only 3D vector
      using VEC2 = ROOT::Math::XYVectorF; // spatial only 2D vector
      using SVEC = ROOT::Math::SVector<double,2>; // 2D algebraic vector
      using SMAT = ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2>>; // 2D symmetric matrix
      TwoDPoint() {}
      TwoDPoint(SVEC const& pos, SMAT const& cov) : pos_(pos), cov_(cov) {} // explicit construct
      TwoDPoint(VEC2 const& pos,UVVar const& uvvar);// construct from a point and resolution information
      TwoDPoint(VEC2 const& pos,VEC2 const& udir, float uvar, float vvar); // construct from a point, semi-major direction cosine to the X axis, and semi-major and minor variances
      TwoDPoint(VEC3 const& pos,VEC3 const& udir, float uvar, float vvar); // construct from a point, semi-major direction, and semi-major and minor variances.  This is used for ComboHit
      // accessors
      auto const& pos() const { return pos_; }
      auto const& cov() const { return cov_; }
      VEC2 pos2() const { return VEC2(pos_[0],pos_[1]); } // physical vector
      VEC3 pos3(float zval=0.0) const { return VEC3(pos_[0],pos_[1],zval); } // physical vector; Z coordinate is meaningless
      // diagonzlied covariance, with angle, used to create new ComboHits
      UVVar uvRes() const;
      void print(std::ostream& os) const;
    private:
      SVEC pos_; // position in XY cartesian coordinates represented as an algebraic vector
      SMAT cov_; // covariance matrix in XY cartesian coordinates
  };
}
std::ostream& operator << (std::ostream& ost, mu2e::TwoDPoint const& pt);
#endif
