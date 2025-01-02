#ifndef Mu2eKinKal_WireHitState_hh
#define Mu2eKinKal_WireHitState_hh
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
#include "Offline/Mu2eKinKal/inc/KKSHFlag.hh"
#include "Offline/Mu2eKinKal/inc/WHSMask.hh"
#include <stdexcept>
#include <iostream>
#include <array>

namespace mu2e {
  // struct describing wire hit internal state
  struct WireHitState {
    enum State { unusable=-3, inactive=-2, left=-1, null=0, right=1};  // drift state
    enum QualityIndex { bkg=0, sign=1, drift=2, chi2=3};// quality info index
    State state_ = null;
    StrawHitUpdaters::algorithm algo_ = StrawHitUpdaters::unknown; // algorithm used to set this state
    bool frozen_ = false; // if set, state not allowed to change during update
    KKSHFlag flag_; // flags for KKStrawHit
    using QType = std::array<double,4>;
    QType quality_ = {-1.0,-1.0,-1.0,-1.0}; // algorithm-dependent, dimensionless quality of this state assignment
// convenience functions
    bool frozen() const { return frozen_; }
    bool wireConstraint() const { return state_ == null; }
    bool driftConstraint() const { return state_ == left || state_ == right; }
    bool isInactive() const { return state_ == inactive; }
    bool active() const { return state_ > inactive; }
    bool usable() const { return state_ > unusable; }
    bool updateable(StrawHitUpdaters::algorithm algo) const { return usable() && ( (!frozen_) || algo_ == algo); } // allow algorithms to update themselves, even if frozen
    auto const& timeConstraint() const { return flag_; }
    bool constrainTOT() const { return flag_.hasAnyProperty(KKSHFlag::tot); }
    bool constrainDriftDt() const { return flag_.hasAnyProperty(KKSHFlag::driftdt); }
    bool constrainAbsDriftDt() const { return flag_.hasAnyProperty(KKSHFlag::absdrift); }
    bool constrainLong() const { return flag_.hasAnyProperty(KKSHFlag::longval); }
    bool nullDriftVar() const { return flag_.hasAnyProperty(KKSHFlag::nhdrift); }
    bool operator == (WireHitState const& whstate) const { return state_ == whstate.state_; }
    bool operator != (WireHitState const& whstate) const { return state_ != whstate.state_; }
    bool operator == (WireHitState::State state) const { return state_ == state; }
    bool operator != (WireHitState::State state) const { return state_ != state; }
    double lrSign() const {
      switch (state_) {
        case left:
          return -1.0;
        case right:
          return 1.0;
        default:
          return 0.0;
      }
    }
    bool isIn(WHSMask const& whsmask) const;
    WireHitState(State state,StrawHitUpdaters::algorithm algo,KKSHFlag const& flag) :
      state_(state), algo_(algo), flag_(flag) {}
    WireHitState(State state = inactive,StrawHitUpdaters::algorithm algo=StrawHitUpdaters::none) :
      state_(state), algo_(algo) {}
    double quality(QualityIndex qid) const { return quality_[static_cast<size_t>(qid)]; }
   // utility functions to convert strings to values
    static std::vector<std::string> StateNames_;
  };
  std::ostream& operator <<(std::ostream& ost, WireHitState const& whs);
}
#endif
