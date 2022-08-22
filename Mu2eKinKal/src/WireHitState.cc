#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
namespace mu2e {
  std::ostream& operator <<(std::ostream& ost, WireHitState const& whs) {
    if(whs.frozen()) ost << " Frozen";
    switch (whs.state_) {
      case WireHitState::unusable:
        ost << " Unusable ";
        break;
      case WireHitState::inactive:
        ost << " Inactive ";
        break;
      case WireHitState::left:
        ost << " Left ";
        break;
      case WireHitState::right:
        ost << " Right ";
        break;
      case WireHitState::null:
        ost << " Null ";
        break;
      default:
         ost << " Unknown ";
        break;
   }
    return ost;
  }
}
