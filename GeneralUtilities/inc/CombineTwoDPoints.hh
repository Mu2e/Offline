//
//  Utility classes for combining 2D points with covariance matrices
//  Original author: Dave Brown (LBNL) 4/7/2023
//
#ifndef GeneralUtilities_CombineTwoDPoints_hh
#define GeneralUtilities_CombineTwoDPoints_hh
#include "Offline/GeneralUtilities/inc/TwoDPoint.hh"
#include "Offline/GeneralUtilities/inc/TwoDWeight.hh"
#include <map>
#include <vector>
#include <iostream>

namespace mu2e {
  class CombineTwoDPoints{
    struct CWT {
      TwoDWeight wt_{}; // weight
      double dchi0_{0}; // chisquared contribution
      CWT() = default;
      CWT(TwoDWeight const& wt, double dchi0) : wt_(wt), dchi0_(dchi0) {}
    };

    public:
      CombineTwoDPoints(float intrinsicvar=0.0) : ivar_(intrinsicvar) {} // empty constructor
      // construct from a vector of points; each points index into the vector is its key
      CombineTwoDPoints(std::vector<TwoDPoint> const& points,float intrinsicvar=0.0);
      // add a point (running average) with a key
      void addPoint(TwoDPoint const& point, size_t key);
      // remove a point by its key
      void removePoint(size_t key);
      // accessors; aggregate point is lazy-evaluated due to inversion
      TwoDPoint const& point() const;
      // return a point by its key
      TwoDPoint point(size_t key) const {return wts_.at(key).wt_.point(ivar_); }
      // aggregate weight
      auto const& weight() const { return wt_; }
      // individual element weight by key
      auto const& weight(size_t key) const {return wts_.at(key).wt_; }
      auto nPoints() const { return wts_.size(); }
      unsigned nDOF() const { return 2*(wts_.size()-1); }
      auto const& weights() const { return wts_; }
      // lazy-evaluated chisquared
      double chisquared() const;
      double consistency() const;
      // compute the chisquared from this object to an external point (not yet added)
      double dChi2(TwoDPoint const& point) const;
      // compute chisquared contribution of a point already included
      double dChi2(size_t key) const;
      void print(std::ostream& os) const;
    private:
      double ivar_ = 0.0; // intrinsic variance
      mutable bool ptcalc_ = false, chicalc_ = false;
      mutable TwoDPoint point_; // current combined point
      TwoDWeight wt_; // current combined weight
      double chi0_ = 0.0;
      mutable double chisq_ = 0.0; // chisquared of current combination
      std::map<size_t,CWT> wts_; // weights used in this combo
  };
}
#endif
