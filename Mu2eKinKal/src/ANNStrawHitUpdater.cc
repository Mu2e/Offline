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
    std::string freeze = std::get<3>(annshuconfig);
    freeze_ = WHSMask(freeze);
    diag_ = std::get<4>(annshuconfig);
    if(diag_ > 0)
      std::cout << "ANNStrawHitUpdater weighs" << std::get<0>(annshuconfig) << " cut " << mvacut_ << " null doca " << nulldoca_ << " freezeing " << freeze_ << std::endl;
    if(diag_ > 1)mva_->showMVA();
  }

  std::string const& ANNStrawHitUpdater::configDescription() {
    static std::string descrip( "Weight file, ANN cut, null hit doca, states to freeze, diag level");
    return descrip;
  }

  WireHitState ANNStrawHitUpdater::wireHitState(WireHitState const& input, ClosestApproachData const& tpdata, DriftInfo const& dinfo, ComboHit const& chit) const {
    WireHitState whstate = input;
    if(input.updateable()){
      whstate.algo_ = StrawHitUpdaters::ANN;
      if(nulldoca_ > 0.0)
        whstate.nulldvar_ = nulldoca_*nulldoca_/3.0; // assumes a flat distribution over [-nulldoca_,nulldoca_]
      else
        // interpret negative nulldoca as the minimum drift distance
        whstate.nulldvar_ = std::max(nulldoca_*nulldoca_,dinfo.driftDistance_*dinfo.driftDistance_)/3.0;

      // invoke the ANN
      std::vector<Float_t> pars(5,0.0);
      // this order is given by the training
      pars[0] = fabs(tpdata.doca());
      pars[1] = dinfo.driftDistance_;
      pars[2] = (dinfo.driftDistance_ -fabs(tpdata.doca()))/sqrt(tpdata.docaVar() + pow(dinfo.driftDistanceError_,2));
      pars[3] = chit.driftTime();
      pars[4] = 1000.0*chit.energyDep();
      float mvaout = mva_->evalMVA(pars);
      whstate.quality_ = mvaout;
      if(mvaout > mvacut_){
        whstate.state_ = tpdata.doca() > 0.0 ? WireHitState::right : WireHitState::left;
      } else {
        whstate.state_ = WireHitState::null;
      }
      whstate.frozen_ = whstate.isIn(freeze_);
    }
    return whstate;
  }
}
