#include "Offline/Mu2eKinKal/inc/CAStrawHitUpdater.hh"
#include <cmath>

namespace mu2e {
  using KinKal::ClosestApproachData;
  using KinKal::VEC3;
  CAStrawHitUpdater::    CAStrawHitUpdater(CASHUConfig const& cashuconfig) {
    mindoca_ = std::get<0>(cashuconfig);
    maxdoca_ = std::get<1>(cashuconfig);
    mindt_ = std::get<2>(cashuconfig);
    maxdt_ = std::get<3>(cashuconfig);
    std::string freeze = std::get<4>(cashuconfig);
    freeze_ = WHSMask(freeze);
    std::cout << "CAStrawHitUpdater doca range [" << mindoca_ << "," << maxdoca_ << "] dt range [" << mindt_ << "," << maxdt_ << "] Freezing " << freeze_ << std::endl;
  }


  WireHitState CAStrawHitUpdater::wireHitState(WireHitState const& input, ClosestApproachData const& tpdata) const {
    WireHitState whstate = input;
    if(input.updateable()){
      whstate.algo_ = StrawHitUpdaters::CA;
      double ddoca = std::min(mindoca_,2.4); // limit to the effective straw radius
      whstate.nulldvar_ = ddoca*ddoca/3.0;
      double absdoca = fabs(tpdata.doca());
      if(tpdata.deltaT() < maxdt_ && tpdata.deltaT() > mindt_) {
        if( absdoca < maxdoca_ && absdoca > mindoca_){  // in the sweet spot: use the DOCA to sign the ambiguity
          whstate.state_ = tpdata.doca() > 0.0 ? WireHitState::right : WireHitState::left;
        } else if(absdoca < mindoca_){
          whstate.state_ = WireHitState::null;
        } else {
          whstate.state_ = WireHitState::inactive;
        }
      } else {
        whstate.state_ = WireHitState::inactive;
      }
      whstate.frozen_ = whstate.isIn(freeze_);
    }
    return whstate;
  }
}
