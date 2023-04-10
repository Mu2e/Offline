//
//  Utility classes for combining 2D points with covariance matrices
//  Original author: Dave Brown (LBNL) 4/7/2023
//
#ifndef GeneralUtilities_CombineTwoDPoints_hh
#define GeneralUtilities_CombineTwoDPoints_hh
#include "Offline/GeneralUtilities/inc/TwoDPoint.hh"
#include <iostream>

namespace mu2e {
  class TwoDWeight{
    public:
      using SVEC = ROOT::Math::SVector<double,2>; // 2D algebraic vector
      using SMAT = ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2>>; // 2D covariance matrix
      TwoDWeight() {}
      TwoDWeight(SVEC const& wtpos, SMAT const& wt) : wtpos_(wtpos), wt_(wt) {} // explicit constructor
      TwoDWeight(TwoDPoint const& pos,float ivar=0.0); // construct from a position
      // accessors
      auto const& weight() const { return wt_; }
      auto const& weightPosition() const { return wtpos_; }
      auto & wt() const { return wt_; }
      auto & wtPos() const { return wtpos_; }

      // combine with another weight
      TwoDWeight& operator +=(TwoDWeight const& other);
      // re-create the original point
      TwoDPoint point() const;
      void print(std::ostream& os) const;
   private:
      SVEC wtpos_; // weighted position
      SMAT wt_; // weight
  };

  class CombinedTwoDPoints{
    public:
      // construct from a point.  Optionally define an intrinsic variance
      CombinedTwoDPoints(TwoDPoint const& point,float intrinsicvar=0.0);
      // construct from a vector of points
      CombinedTwoDPoints(std::vector<TwoDPoint> points,float intrinsicvar=0.0);
      // add a point (running average)
      void addPoint(TwoDPoint const& point);
      // allow removing an existing point: not sure how to identify points TODO
      // accessors; point is lazy-evaluated
      TwoDPoint const& point() const;
      auto const& weight() const { return wt_; }
      auto nPoints() const { return wts_.size(); }
      unsigned nDOF() const { return 2*(wts_.size()-1); }
      // lazy-evaluated chisquared
      double consistency() const;
      double chisquared() const;
      // compute the chisquared from this object to a point
      double dChi2(TwoDPoint const& point) const;
      void print(std::ostream& os) const;
    private:
      double ivar_ = 0.0; // intrinsic variance
      mutable bool ptcalc_, chicalc_;
      mutable TwoDPoint point_; // current point
      TwoDWeight wt_; // current weight sum
      double chi0_ = 0.0;
      mutable double chisq_; // chisquared of current combination
      std::vector<TwoDWeight> wts_; // weights used in this combination
  };
}
std::ostream& operator << (std::ostream& ost, mu2e::TwoDWeight const& pt);
#endif
