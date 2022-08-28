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
    mvacut_ = std::get<1>(annshuconfig);
    nulldoca_ = std::get<2>(annshuconfig);
    std::string freeze = std::get<3>(annshuconfig);
    freeze_ = WHSMask(freeze);
    std::cout << "ANNStrawHitUpdater " << " anncut " << mvacut_ << " null doca" << nulldoca_ << " freezeing " << freeze_ << std::endl;
    mva_->initMVA();
    mva_->showMVA();
  }

  WireHitState ANNStrawHitUpdater::wireHitState(WireHitState const& input, ClosestApproachData const& tpdata, DriftInfo const& dinfo, ComboHit const& chit) const {
    WireHitState whstate = input;
    if(input.updateable()){
      std::vector<Float_t> pars(8,0.0);
      // this order is given by the training
      pars[0] = fabs(tpdata.doca());
      pars[1] = dinfo.driftDistance_;
      pars[2] = chit.driftTime();
      pars[3] = 1000.0*chit.energyDep();
      pars[4] = tpdata.docaVar();
      pars[5] = fabs(tpdata.dirDot());
      pars[6] = fabs(tpdata.particlePoca().Vect().Z());
      pars[7] = tpdata.particlePoca().Vect().Rho();
      float mvaout = mva_->evalMVA(pars);
      whstate.algo_ = StrawHitUpdaters::ANN;
      whstate.nulldoca_ = nulldoca_;
      if(mvaout > mvacut_){
        whstate.state_ = tpdata.doca() > 0.0 ? WireHitState::right : WireHitState::left;
      } else {
        whstate.state_ = WireHitState::null;
      } // add MVA for inactivating hits TODO
      whstate.frozen_ = whstate.isIn(freeze_);
    }
    return whstate;
  }
}
