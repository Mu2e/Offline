#include "Offline/Mu2eKinKal/inc/PanelDiagSHU.hh"
#include <cmath>

namespace mu2e {
  PanelDiagSHU::PanelDiagSHU(Config const& config) {
    mask_ = StrawIdMask(std::get<0>(config));
    sid_ = StrawId(std::get<1>(config));
    diag_ = std::get<2>(config);
    if(diag_ > 0)std::cout << "PanelDiagSHU " << mask_.levelName() << " " << sid_ << std::endl;
  }

  WireHitState PanelDiagSHU::wireHitState(WireHitState const& input, StrawId const& id) const {
    WireHitState whstate = input;
    if(input.updateable(StrawHitUpdaters::PanelDiag)){
      if (mask_.equal(id,sid_)){
        whstate.state_ = WireHitState::inactive;
        whstate.algo_ = StrawHitUpdaters::PanelDiag;
        whstate.frozen_ = true;
        if (diag_ > 1)std::cout << "PanelDiagSHU set hit " << id << " inactive" << std::endl;
      }
    }
    return whstate;
  }

  std::string const& PanelDiagSHU::configDescription() {
    static std::string descrip( "StrawId mask, StrawId to set inactive, diag level");
    return descrip;
  }

}
