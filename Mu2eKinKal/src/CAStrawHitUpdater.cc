#include "Offline/Mu2eKinKal/inc/CAStrawHitUpdater.hh"
#include <cmath>

namespace mu2e {
  using KinKal::ClosestApproachData;
  using KinKal::VEC3;
  CAStrawHitUpdater::    CAStrawHitUpdater(CASHUConfig const& cashuconfig) {
    mindoca_ = std::get<0>(cashuconfig);
    maxdoca_ = std::get<1>(cashuconfig);
    minrdrift_ = std::get<2>(cashuconfig);
    maxrdrift_ = std::get<3>(cashuconfig);
    std::string freeze = std::get<4>(cashuconfig);
    freeze_ = WHSMask(freeze);
    std::cout << "CAStrawHitUpdater doca range [" << mindoca_ << "," << maxdoca_ << "] rdrift range [" << minrdrift_ << "," << maxrdrift_ << "] Freezing " << freeze_ << std::endl;
  }

  WireHitState CAStrawHitUpdater::wireHitState(WireHitState const& input, ClosestApproachData const& tpdata,DriftInfo const& dinfo) const {
    WireHitState whstate = input;
    if(input.updateable()){
      whstate.algo_ = StrawHitUpdaters::CA;
      double rd = std::min(minrdrift_,2.4); // limit to the effective straw radius
      whstate.nulldvar_ = rd*rd/3.0;
      double absdoca = fabs(tpdata.doca());
      if(dinfo.driftDistance_ < maxrdrift_ && dinfo.driftDistance_ > minrdrift_ &&
          absdoca < maxdoca_ && absdoca > mindoca_){  // in the sweet spot: use the DOCA to sign the ambiguity
        whstate.state_ = tpdata.doca() > 0.0 ? WireHitState::right : WireHitState::left;
      } else if(dinfo.driftDistance_ < minrdrift_){
        whstate.state_ = WireHitState::null;
      } else {
        whstate.state_ = WireHitState::inactive;
      }
      whstate.frozen_ = whstate.isIn(freeze_);
    }
    return whstate;
  }
}
