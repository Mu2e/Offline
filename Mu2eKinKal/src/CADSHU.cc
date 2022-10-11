#include "Offline/Mu2eKinKal/inc/CADSHU.hh"
#include <cmath>

namespace mu2e {
  using KinKal::ClosestApproachData;
  using KinKal::VEC3;
  CADSHU::CADSHU(Config const& config) {
    maxdoca_ = std::get<0>(config);
    double maxdocaerr = std::get<1>(config);
    maxdvar_ = maxdocaerr*maxdocaerr;
    minrdrift_ = std::get<2>(config);
    maxrdrift_ = std::get<3>(config);
    std::string allowed = std::get<4>(config);
    allowed_ = WHSMask(allowed);
    std::string freeze = std::get<5>(config);
    freeze_ = WHSMask(freeze);
    // set the null hit variance
    double rd = std::min(minrdrift_,2.4); // limit to the effective straw radius. this value should be configurable TODO
    nulldvar_ = rd*rd/3.0;
    diag_ = std::get<6>(config);
    if(diag_ > 0)std::cout << "CADSHU max doca, doca error " << maxdoca_ << " " << maxdocaerr
      << " rdrift range [" << minrdrift_ << "," << maxrdrift_ << "] Allowing " << allowed_ << " Freezing " << freeze_ << std::endl;
  }

  WireHitState CADSHU::wireHitState(WireHitState const& input, ClosestApproachData const& tpdata,DriftInfo const& dinfo) const {
    WireHitState whstate = input;
    if(input.updateable(StrawHitUpdaters::CAD)){
      double absdoca = fabs(tpdata.doca());
      if(dinfo.driftDistance_ < maxrdrift_ && absdoca < maxdoca_ && tpdata.docaVar() > 0.0 && tpdata.docaVar() < maxdvar_ ){
        if(dinfo.driftDistance_ > minrdrift_){
          // in the sweet spot: use the DOCA to sign the ambiguity
          if(allowed_.hasAnyProperty(WHSMask::drift)) {
            whstate.state_ = tpdata.doca() > 0.0 ? WireHitState::right : WireHitState::left;
            whstate.algo_ = StrawHitUpdaters::CAD;
          }
        } else if(allowed_.hasAnyProperty(WHSMask::null)) {
          whstate.state_ = WireHitState::null;
          whstate.algo_ = StrawHitUpdaters::CAD;
          whstate.nulldvar_ = nulldvar_;
        }
      } else if(allowed_.hasAnyProperty(WHSMask::inactive)) {
        whstate.state_ = WireHitState::inactive;
        whstate.algo_ = StrawHitUpdaters::CAD;
      }
      if(whstate.algo_ == StrawHitUpdaters::CAD)whstate.frozen_ = whstate.isIn(freeze_);
    }
    return whstate;
  }

  std::string const& CADSHU::configDescription() {
    static std::string descrip( "Maximum DOCA to use hit, Maximum DOCA error to use hit, Minimum rdrift to set LR ambiguity, Maximum rdrift to use hit, allowed states, States to freeze, diag level");
    return descrip;
  }

}
