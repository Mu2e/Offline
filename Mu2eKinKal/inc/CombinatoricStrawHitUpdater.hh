#ifndef Mu2eKinKal_CombinatoricStrawHitUpdater_hh
#define Mu2eKinKal_CombinatoricStrawHitUpdater_hh
//
//  StrawHitGroup updating using an exhaustive combinatoric algorithm, following the BTrk PanelAmbigResolver algorithm
//
#include "KinKal/Detector/WireHitStructs.hh"

#include <tuple>
#include <vector>

namespace mu2e {
  using KinKal::WireHitState;
  using WHSCOL = std::vector<WireHitState>;
  using CHI2WHS = std::tuple<double,WHSCOL>;
  using CHI2WHSCOL = std::vector<CHI2WHS>;

  class CombinatoricStrawHitUpdater {
    public:
      struct CHI2Comp { // comparator to sort hit states by chisquared value
        bool operator()(CHI2WHS const& a, CHI2WHS const& b)  const {
          return std::get<0>(a) < std::get<0>(b);
        }
      };
      CombinatoricStrawHitUpdater(double inactivep, double nullp, double mindchi2) :
        inactivep_(inactivep), nullp_(nullp), mindchi2_(mindchi2),
        allowed_{WireHitState::inactive, WireHitState::left, WireHitState::null, WireHitState::right} {}
      WHSCOL selectBest(CHI2WHSCOL& chi2s) const; // find the best configuration given the total chisq for each
      double penalty(WireHitState const& whs) const; // compute the penalty for each hit in a given state
      double inactivePenalty() const { return inactivep_;}
      double nullPenalty() const { return nullp_;}
      auto const& allowed() const { return allowed_; }
      double minDeltaChi2() const { return mindchi2_; }
    private:
      double inactivep_; // chisquared penalty for inactive hits
      double nullp_; // chisquared penalty for null hits
      double mindchi2_; // minimum chisquared separation to consider 'significant'
      WHSCOL allowed_; // allowed states
  };
}
#endif
