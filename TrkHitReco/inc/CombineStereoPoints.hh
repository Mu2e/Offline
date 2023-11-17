//
//  Utility classes for combining StereoPoints.  This is similar to Combine2DPoints, but allows the
//  stereo hits to have an internal direction (du/dw, dv/dw)
//  Original author: Dave Brown (LBNL) 11/15/2023
//
#ifndef TrkHitReco_CombineStereoPoints_hh
#define TrkHitReco_CombineStereoPoints_hh
#include "Offline/TrkHitReco/inc/StereoPoint.hh"
#include "Offline/GeneralUtilities/inc/TwoDWeight.hh"
#include <map>
#include <vector>
#include <iostream>

namespace mu2e {
  class CombineStereoPoints{
    struct CWT {
      TwoDWeight wt_; // weight
      double z_; // z position of this hit
      double dchi0_; // chisquared contribution
      CWT(TwoDWeight const& wt, double z, double dchi0) : wt_(wt), z_(z), dchi0_(dchi0) {}
    };

    public:
      CombineStereoPoints(float intrinsicvar=0.0) : ivar_(intrinsicvar) {} // empty constructor
      // add a point (running average) with a key
      void addPoint(StereoPoint const& point, size_t key);
      // remove a point by its key
      void removePoint(size_t key);
      // return a point by its key
      StereoPoint point(size_t key) const;
      // return aggreagate point
      StereoPoint const& point() const;
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
      double dChi2(StereoPoint const& point) const;
      // compute chisquared contribution of a point already included
      double dChi2(size_t key) const;
      void print(std::ostream& os) const;
    private:
      double ivar_ = 0.0; // intrinsic variance
      mutable bool ptcalc_ = false, chicalc_ = false;
      TwoDWeight wt_; // current combined weight
      mutable StereoPoint point_; // current combined point (lazy evaluated)
      double chi0_ = 0.0;
      mutable double chisq_ = 0.0; // chisquared of current combination
      std::map<size_t,CWT> wts_; // weights used in this combo
  };
}
#endif
