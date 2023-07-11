//
//  weight-space representation of a 2D point constraint
//
#ifndef GeneralUtilities_TwoDWeight_hh
#define GeneralUtilities_TwoDWeight_hh
#include "Offline/GeneralUtilities/inc/TwoDPoint.hh"
namespace mu2e{
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
      TwoDWeight& operator -=(TwoDWeight const& other);
      // re-create the original point
      TwoDPoint point(float ivar=0.0) const;
      void print(std::ostream& os) const;
    private:
      SVEC wtpos_; // weighted position
      SMAT wt_; // weight
  };
}
std::ostream& operator << (std::ostream& ost, mu2e::TwoDWeight const& pt);
#endif
