#ifndef Mu2eKinKal_WireHitState_hh
#define Mu2eKinKal_WireHitState_hh
#include <stdexcept>
#include <iostream>

namespace mu2e {
  // struct describing wire hit internal state
  struct WireHitState {
    enum State { forcedinactive=-3, inactive=-2, left=-1, null=0, right=1};  // state description
    State state_;
    bool useDrift() const { return state_ == left || state_ == right; }
    bool active() const { return state_ > inactive; }
    bool usable() const { return state_ > forcedinactive; }
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
    WireHitState(State state = inactive) : state_(state) {}
    WireHitState& operator = (State state) { state_ = state; return *this; }
  };
  std::ostream& operator <<(std::ostream& ost, WireHitState const& whs);
}
#endif
