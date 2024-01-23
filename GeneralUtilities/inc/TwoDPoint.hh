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
  class TwoDPoint {
    public:
      using VEC3 = ROOT::Math::XYZVectorF; // spatial only 3D vector
      using VEC2 = ROOT::Math::XYVectorF; // spatial only 2D vector
      using SVEC = ROOT::Math::SVector<double,2>; // 2D algebraic vector
      using SMAT = ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2>>; // 2D symmetric matrix
      TwoDPoint() {}
      TwoDPoint(SVEC const& pos, SMAT const& cov); // explicit construct
      TwoDPoint(VEC2 const& pos,VEC2 const& udir, float uvar, float vvar); // construct from a point, semi-major direction cosine to the X axis, and semi-major and minor variances
      TwoDPoint(VEC3 const& pos,VEC3 const& udir, float uvar, float vvar); // construct from a point, semi-major direction, and semi-major and minor variances.  This is used for ComboHit
      // accessors
      auto const& pos() const { return pos_; }
      auto const& cov() const { return cov_; }
      VEC2 pos2() const { return VEC2(pos_[0],pos_[1]); } // physical vector
      VEC3 pos3(float zval=0.0) const { return VEC3(pos_[0],pos_[1],zval); } // physical vector; Z coordinate is meaningless
      // eigenvalues and eigenvectors of the covariance, used to create new ComboHits
      VEC2 const& udir() const { return udir_;}
      VEC2 vdir() const { return VEC2(-udir_.Y(),udir_.X()); }
      float uvar() const;
      float vvar() const;
      void print(std::ostream& os) const;
    private:
      SVEC pos_; // position in XY cartesian coordinates represented as an algebraic vector
      SMAT cov_; // covariance matrix in XY cartesian coordinates
      VEC2 udir_; // direction of covariance semi-major axis
  };
}
std::ostream& operator << (std::ostream& ost, mu2e::TwoDPoint const& pt);
#endif
