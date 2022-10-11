#include "Offline/Mu2eKinKal/inc/BkgStrawHitUpdater.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include <cmath>

namespace mu2e {
  using KinKal::ClosestApproachData;
  using KinKal::VEC3;
  BkgStrawHitUpdater::BkgStrawHitUpdater(BkgSHUConfig const& bkgshuconfig) {
    mva_  = new MVATools(std::get<0>(bkgshuconfig));
    mva_->initMVA();
    mvacut_ = std::get<1>(bkgshuconfig);
    std::string freeze = std::get<2>(bkgshuconfig);
    diag_ = std::get<3>(bkgshuconfig);
    freeze_ = WHSMask(freeze);
    if(diag_ > 0)std::cout << "BkgStrawHitUpdater weights " << std::get<0>(bkgshuconfig) << " cut " << mvacut_ << " freeze " << freeze_ << std::endl;
    if(diag_ > 1)mva_->showMVA();
  }

  WireHitState BkgStrawHitUpdater::wireHitState(WireHitState const& input, ClosestApproachData const& tpdata, DriftInfo const& dinfo, ComboHit const& chit) const {
    WireHitState whstate = input;
    if(input.updateable(StrawHitUpdaters::Bkg)){
      std::vector<Float_t> pars(10,0.0);
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
      // compare the delta-t based U position with the fit U position; requires relative end
      double endsign = 2.0*(chit.driftEnd()-0.5);
      double upos = -endsign*tpdata.sensorDirection().Dot(tpdata.sensorPoca().Vect() - chit.centerPos());
      pars[7] = fabs(chit.wireDist()-upos)/chit.wireRes();
      pars[8] = chit.wireRes();
      pars[9] = tpdata.particlePoca().Vect().Rho();
      float mvaout = mva_->evalMVA(pars);
      if(diag_ > 2){
        whstate.algo_  = StrawHitUpdaters::Bkg;
        whstate.quality_ = mvaout;
      }
      if(mvaout < mvacut_){
        whstate.algo_  = StrawHitUpdaters::Bkg;
        whstate.state_ = WireHitState::inactive;
        whstate.frozen_ = whstate.isIn(freeze_);
        whstate.quality_ = mvaout;
      }
    }
    return whstate;
  }
  std::string const& BkgStrawHitUpdater::configDescription() {
    static std::string descrip( "Weight file, ANN cut, states to freeze, diag level");
    return descrip;
  }
}
