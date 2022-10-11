#include "Offline/Mu2eKinKal/inc/CAStrawHitUpdater.hh"
#include <cmath>

namespace mu2e {
  using KinKal::ClosestApproachData;
  using KinKal::VEC3;
  CAStrawHitUpdater::    CAStrawHitUpdater(CASHUConfig const& cashuconfig) {
    maxdoca_ = std::get<0>(cashuconfig);
    minrdrift_ = std::get<1>(cashuconfig);
    maxrdrift_ = std::get<2>(cashuconfig);
    std::string allowed = std::get<3>(cashuconfig);
    allowed_ = WHSMask(allowed);
    std::string freeze = std::get<4>(cashuconfig);
    freeze_ = WHSMask(freeze);
    // set the null hit variance
    double rd = std::min(minrdrift_,2.4); // limit to the effective straw radius. this value should be configurable TODO
    nulldvar_ = rd*rd/3.0;
    diag_ = std::get<5>(cashuconfig);
    if(diag_ > 0)std::cout << "CAStrawHitUpdater max doca " << maxdoca_ << " rdrift range [" << minrdrift_ << "," << maxrdrift_ << "] Allowing " << allowed_ << " Freezing " << freeze_ << std::endl;
  }

  WireHitState CAStrawHitUpdater::wireHitState(WireHitState const& input, ClosestApproachData const& tpdata,DriftInfo const& dinfo) const {
    WireHitState whstate = input;
    if(input.updateable(StrawHitUpdaters::CA)){
      double absdoca = fabs(tpdata.doca());
      if(dinfo.driftDistance_ < maxrdrift_ && absdoca < maxdoca_){
        if(dinfo.driftDistance_ > minrdrift_){
          // in the sweet spot: use the DOCA to sign the ambiguity
          if(allowed_.hasAnyProperty(WHSMask::drift)) {
            whstate.state_ = tpdata.doca() > 0.0 ? WireHitState::right : WireHitState::left;
            whstate.algo_ = StrawHitUpdaters::CA;
          }
        } else if(allowed_.hasAnyProperty(WHSMask::null)) {
          whstate.state_ = WireHitState::null;
          whstate.algo_ = StrawHitUpdaters::CA;
          whstate.nulldvar_ = nulldvar_;
        }
      } else if(allowed_.hasAnyProperty(WHSMask::inactive)) {
        whstate.state_ = WireHitState::inactive;
        whstate.algo_ = StrawHitUpdaters::CA;
      }
      if(whstate.algo_ == StrawHitUpdaters::CA)whstate.frozen_ = whstate.isIn(freeze_);
    }
    return whstate;
  }

  std::string const& CAStrawHitUpdater::configDescription() {
    static std::string descrip( "Maximum DOCA to use hit, Minimum rdrift to set LR ambiguity, Maximum rdrift to use hit, allowed states, States to freeze, diag level");
    return descrip;
  }

}
