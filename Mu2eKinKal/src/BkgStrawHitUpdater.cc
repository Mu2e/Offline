#include "Offline/Mu2eKinKal/inc/BkgStrawHitUpdater.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include <cmath>

namespace mu2e {
  using KinKal::ClosestApproachData;
  using KinKal::VEC3;
  BkgStrawHitUpdater::BkgStrawHitUpdater(BkgSHUConfig const& bkgshuconfig) {
    mva_  = new MVATools(std::get<0>(bkgshuconfig));
    mvacut_ = std::get<1>(bkgshuconfig);
    std::string freeze = std::get<2>(bkgshuconfig);
    freeze_ = WHSMask(freeze);
    std::cout << "BkgStrawHitUpdater " << " bkgcut " << mvacut_ << " freeze " << freeze_ << std::endl;
    mva_->initMVA();
    mva_->showMVA();
  }

  WireHitState BkgStrawHitUpdater::wireHitState(WireHitState const& input, ClosestApproachData const& tpdata, DriftInfo const& dinfo, ComboHit const& chit) const {
    WireHitState whstate = input;
    if(input.updateable()){
      std::vector<Float_t> pars(11,0.0);
      // this order is given by the training
      pars[0] = fabs(tpdata.doca());
      pars[1] = dinfo.driftDistance_;
      pars[2] = chit.driftTime();
      pars[3] = 1000.0*chit.energyDep();
      pars[4] = tpdata.docaVar();
      pars[5] = fabs(tpdata.dirDot());
      pars[6] = fabs(tpdata.particlePoca().Vect().Z());
      pars[7] = tpdata.particlePoca().Vect().Rho();
      // compare the delta-t based U position with the fit U position
      double endsign = 2.0*(chit.driftEnd()-0.5);
      double upos = -endsign*tpdata.sensorDirection().Dot(tpdata.sensorPoca().Vect() - chit.centerPos());
      pars[8] = fabs(chit.wireDist()-upos)/chit.wireRes();
      pars[9] = chit.wireRes();
      pars[10] = dinfo.LorentzAngle_;
      float mvaout = mva_->evalMVA(pars);
      if(mvaout < mvacut_){
        whstate.algo_  = StrawHitUpdaters::Bkg;
        whstate.state_ = WireHitState::inactive;
        whstate.frozen_ = whstate.isIn(freeze_);
      }
    }
    return whstate;
  }
}
