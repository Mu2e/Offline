#include "Offline/Mu2eKinKal/inc/ANNStrawHitUpdater.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include <cmath>

namespace mu2e {
  using KinKal::ClosestApproachData;
  using KinKal::VEC3;
  ANNStrawHitUpdater::ANNStrawHitUpdater(ANNSHUConfig const& annshuconfig) {
    mva_  = new MVATools(std::get<0>(annshuconfig));
    mva_->initMVA();
    mvacut_ = std::get<1>(annshuconfig);
    nulldoca_ = std::get<2>(annshuconfig);
    std::string allowed = std::get<3>(annshuconfig);
    allowed_ = WHSMask(allowed);
    std::string freeze = std::get<4>(annshuconfig);
    freeze_ = WHSMask(freeze);
    diag_ = std::get<5>(annshuconfig);
    if(diag_ > 0)
      std::cout << "ANNStrawHitUpdater weighs" << std::get<0>(annshuconfig) << " cut " << mvacut_ << " null doca " << nulldoca_
        << " allowing " << allowed_ << " freezing " << freeze_ << std::endl;
    if(diag_ > 1)mva_->showMVA();
  }

  std::string const& ANNStrawHitUpdater::configDescription() {
    static std::string descrip( "Weight file, ANN cut, null hit doca, allowed states, states to freeze, diag level");
    return descrip;
  }

  WireHitState ANNStrawHitUpdater::wireHitState(WireHitState const& input, ClosestApproachData const& tpdata, DriftInfo const& dinfo, ComboHit const& chit) const {
    WireHitState whstate = input;
    if(input.updateable(StrawHitUpdaters::ANN)){
      // invoke the ANN
      std::vector<Float_t> pars(7,0.0);
      // this order is given by the training
      double derr = sqrt(std::max(tpdata.docaVar(),0.0) + dinfo.driftDistanceError_*dinfo.driftDistanceError_);
      double dvar = derr*derr;
      double totvar = dvar + tpdata.docaVar();
      double udoca = fabs(tpdata.doca());
      pars[0] = udoca;
      pars[1] = dinfo.driftDistance_;
      pars[2] = (dinfo.driftDistance_ - udoca)/sqrt(totvar);
      pars[3] = chit.driftTime();
      pars[4] = 1000.0*chit.energyDep();
      pars[5] = derr;
      pars[6] =4.0*dinfo.driftDistance_*udoca/totvar;
      float mvaout = mva_->evalMVA(pars);
      whstate.quality_ = mvaout;
      if(mvaout > mvacut_){
        if(allowed_.hasAnyProperty(WHSMask::drift)){
          whstate.state_ = tpdata.doca() > 0.0 ? WireHitState::right : WireHitState::left;
          whstate.algo_ = StrawHitUpdaters::ANN;
        }
      } else {
        if(allowed_.hasAnyProperty(WHSMask::null)){
          whstate.state_ = WireHitState::null;
          whstate.algo_ = StrawHitUpdaters::ANN;
          if(nulldoca_ > 0.0)
            whstate.nulldvar_ = nulldoca_*nulldoca_/3.0; // assumes a flat distribution over [-nulldoca_,nulldoca_]
          else
            // interpret negative nulldoca as the minimum drift distance
            whstate.nulldvar_ = std::max(nulldoca_*nulldoca_,dinfo.driftDistance_*dinfo.driftDistance_)/3.0;
        }
      }
      if(whstate.algo_ == StrawHitUpdaters::ANN)whstate.frozen_ = whstate.isIn(freeze_);
    }
    return whstate;
  }
}
