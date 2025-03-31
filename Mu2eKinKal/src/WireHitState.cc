#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "cetlib_except/exception.h"
namespace mu2e {
  std::vector<std::string> WireHitState::StateNames_ = {"Unusable", "Inactive", "Left", "Null","Right"};
  std::ostream& operator <<(std::ostream& ost, WireHitState const& whs) {
    ost << "WireHitState ";
    if(whs.frozen()) ost << " Frozen ";
    ost << WireHitState::StateNames_[whs.state_+3];
    ost << " Flags "  << whs.flag_;
    ost << " algo " << StrawHitUpdaters::name(whs.algo_);
    return ost;
  }

  bool WireHitState::isIn(WHSMask const& whsmask) const {
    switch (state_) {
      case WireHitState::inactive:
        return whsmask.hasAnyProperty(WHSMask::inactive);
      case WireHitState::left: case WireHitState::right:
        return whsmask.hasAnyProperty(WHSMask::drift);
      case WireHitState::null:
        return whsmask.hasAnyProperty(WHSMask::null);
      default:
        return false;
    }
  }
}
