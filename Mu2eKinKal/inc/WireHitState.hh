#ifndef Mu2eKinKal_WireHitState_hh
#define Mu2eKinKal_WireHitState_hh
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
#include "Offline/Mu2eKinKal/inc/WHSMask.hh"
#include <stdexcept>
#include <iostream>

namespace mu2e {
  // struct describing wire hit internal state
  struct WireHitState {
    enum State { unusable=-3, inactive=-2, left=-1, null=0, right=1, leftdt=2,rightdt=3};  // drift state
    enum TOTUse { unused=0, nullonly=1, driftonly=2, all=3}; // how to use TOT time constraint
    enum NullDistVar { rstraw=0, rdrift=1}; // how to assign null variance
    enum QualityIndex { bkg=0, sign=1, drift=2, chi2=3};// quality info index
    State state_ = null;
    StrawHitUpdaters::algorithm algo_ = StrawHitUpdaters::unknown; // algorithm used to set this state
    NullDistVar nulldvar_ = rstraw; // distance variance for null hits
    bool frozen_ = false; // if set, state not allowed to change during update
    TOTUse totuse_ = all;
    std::array<double,4> quality_ = {-1.0,-1.0,-1.0,-1.0}; // algorithm-dependent, dimensionless quality of this state assignment
// convenience functions
    bool frozen() const { return frozen_; }
    bool wireConstraint() const { return state_ == null; }
    bool driftConstraint() const { return state_ == left || state_ == right || state_ == leftdt || state_ == rightdt; }
    bool driftDTConstraint() const { return state_ == leftdt || state_ == rightdt; }
    bool isInactive() const { return state_ == inactive; }
    bool active() const { return state_ > inactive; }
    bool usable() const { return state_ > unusable; }
    bool updateable(StrawHitUpdaters::algorithm algo) const { return usable() && ( (!frozen_) || algo_ == algo); } // allow algorithms to update themselves, even if frozen
    bool constrainTOT() const {
      return totuse_ == all ||
        (wireConstraint() && totuse_ == nullonly) ||
        (driftConstraint() && totuse_ == WireHitState::driftonly);
    }
    bool operator == (WireHitState const& whstate) const { return state_ == whstate.state_; }
    bool operator != (WireHitState const& whstate) const { return state_ != whstate.state_; }
    bool operator == (WireHitState::State state) const { return state_ == state; }
    bool operator != (WireHitState::State state) const { return state_ != state; }
    double lrSign() const {
      switch (state_) {
        case left: case leftdt:
          return -1.0;
        case right: case rightdt:
          return 1.0;
        default:
          return 0.0;
      }
    }
    bool isIn(WHSMask const& whsmask) const;
    WireHitState(State state = inactive,StrawHitUpdaters::algorithm algo=StrawHitUpdaters::none,
        NullDistVar nulldvar=rstraw) : state_(state), algo_(algo) , nulldvar_(nulldvar){}
    // utility functions to convert strings to values
    static TOTUse totUse(std::string const& totuse);
    static NullDistVar nullDistVar(std::string const& ndvar);
    static std::vector<std::string> StateNames_;
    static std::vector<std::string> TOTUseNames_;
    static std::vector<std::string> NullDistVarNames_;
    double quality(QualityIndex qid) const { return quality_[static_cast<size_t>(qid)]; }
  };
  std::ostream& operator <<(std::ostream& ost, WireHitState const& whs);
}
#endif
