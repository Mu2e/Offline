#ifndef Mu2eKinKal_WireHitState_hh
#define Mu2eKinKal_WireHitState_hh
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
#include "Offline/Mu2eKinKal/inc/WHSMask.hh"
#include <stdexcept>
#include <iostream>

namespace mu2e {
  // struct describing wire hit internal state
  struct WireHitState {
    enum State { unusable=-3, inactive=-2, left=-1, null=0, right=1};  // drift state
    enum TOTUse { unused=0, nullonly=1, driftonly=2, all=3}; // how to use TOT time constraint
    State state_;
    StrawHitUpdaters::algorithm algo_; // algorithm used to set this state
    double nulldvar_; // distance variance for null hits
    bool frozen_ = false; // if set, state not allowed to change during update
    TOTUse totuse_ = all;
    double quality_ = -1.0; // algorithm-dependent, dimensionless quality of this state assignment
// convenience functions
    bool frozen() const { return frozen_; }
    bool wireConstraint() const { return state_ == null; }
    bool driftConstraint() const { return state_ == left || state_ == right; }
    bool isInactive() const { return state_ == inactive; }
    bool active() const { return state_ > inactive; }
    bool usable() const { return state_ > unusable; }
    bool updateable(StrawHitUpdaters::algorithm algo) const { return usable() && (!frozen_ || algo_ == algo); } // allow algorithms to update themselves, even if frozen
    double nullDistanceVariance() const { return nulldvar_; }
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
    WireHitState(State state = inactive,StrawHitUpdaters::algorithm algo=StrawHitUpdaters::none,double nulldvar=1.92) : state_(state), algo_(algo), nulldvar_(nulldvar) {}
  };
  std::ostream& operator <<(std::ostream& ost, WireHitState const& whs);
}
#endif
