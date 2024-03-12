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
    std::string flag = std::get<4>(config);
    flag_ = KKSHFlag(flag);
    std::string allowed = std::get<5>(config);
    allowed_ = WHSMask(allowed);
    std::string freeze = std::get<6>(config);
    freeze_ = WHSMask(freeze);
    diag_ = std::get<7>(config);
    if(diag_ > 0)std::cout << "CADSHU max doca, doca error " << maxdoca_ << " " << maxdocaerr
      << " rdrift range [" << minrdrift_ << "," << maxrdrift_ << "] Allowing "
        << allowed_ << " Freezing " << freeze_ << " Flags " << flag
        << std::endl;
  }

  WireHitState CADSHU::wireHitState(WireHitState const& input, ClosestApproachData const& tpdata,DriftInfo const& dinfo) const {
    WireHitState whstate = input;
    bool updated(false);
    if(input.updateable(StrawHitUpdaters::CAD)){
      double absdoca = fabs(tpdata.doca());
      if(dinfo.rDrift_ < maxrdrift_ && absdoca < maxdoca_ && tpdata.docaVar() > 0.0 && tpdata.docaVar() < maxdvar_ ){
        if(dinfo.rDrift_ > minrdrift_){
          // in the sweet spot: use the DOCA to sign the ambiguity
          if(allowed_.hasAnyProperty(WHSMask::drift)) {
            whstate.state_ = tpdata.doca() > 0.0 ? WireHitState::right : WireHitState::left;
            updated = true;
          }
        } else if(allowed_.hasAnyProperty(WHSMask::null)) {
          whstate.state_ = WireHitState::null;
          updated = true;
        }
      } else if(allowed_.hasAnyProperty(WHSMask::inactive)) {
        whstate.state_ = WireHitState::inactive;
        updated = true;
      }
      if(updated){
        whstate.algo_ = StrawHitUpdaters::CAD;
        whstate.flag_ = flag_;
        whstate.frozen_ = whstate.isIn(freeze_);
      }
      if (diag_ > 1)std::cout << "CADSHU set hit " << whstate << std::endl;
    } else if (diag_ > 1) {
      std::cout << "CADSHU skipping hit " << whstate << std::endl;    }
    return whstate;
  }

  std::string const& CADSHU::configDescription() {
    static std::string descrip( "Maximum DOCA to use hit, Maximum DOCA error to use hit, Minimum rdrift to set LR ambiguity, Maximum rdrift to use hit, allowed states, States to freeze, diag level");
    return descrip;
  }

}
