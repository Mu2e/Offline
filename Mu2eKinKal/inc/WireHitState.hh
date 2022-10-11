#ifndef Mu2eKinKal_WireHitState_hh
#define Mu2eKinKal_WireHitState_hh
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
#include "Offline/Mu2eKinKal/inc/WHSMask.hh"
#include <stdexcept>
#include <iostream>

namespace mu2e {
  // struct describing wire hit internal state
  struct WireHitState {
    enum State { unusable=-3, inactive=-2, left=-1, null=0, right=1};  // state description; split null into baddrift and digital
    State state_;
    StrawHitUpdaters::algorithm algo_; // algorithm used to set this state
    double nulldvar_; // distance variance for null hits
    bool frozen_; // if set, state not allowed to change during update
    double quality_; // algorithm-dependent, dimensionless quality of this state assignment
// convenience functions
    bool frozen() const { return frozen_; }
    bool useDrift() const { return state_ == left || state_ == right; }
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
    WireHitState(State state = inactive,StrawHitUpdaters::algorithm algo=StrawHitUpdaters::none,double nulldvar=1.92) : state_(state), algo_(algo), nulldvar_(nulldvar), frozen_(false), quality_(-1.0) {}
  };
  std::ostream& operator <<(std::ostream& ost, WireHitState const& whs);
}
#endif
