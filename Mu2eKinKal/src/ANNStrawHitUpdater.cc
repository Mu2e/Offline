#include "Offline/Mu2eKinKal/inc/ANNStrawHitUpdater.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include <cmath>

namespace mu2e {
  using KinKal::ClosestApproachData;
  using KinKal::VEC3;
  ANNStrawHitUpdater::ANNStrawHitUpdater(ANNSHUConfig const& annshutuple) {
    mva_  = new MVATools(std::get<0>(annshutuple));
    mvacut_ = std::get<1>(annshutuple);
    nulldoca_ = std::get<2>(annshutuple);
    std::string freeze = std::get<3>(annshutuple);
    freeze_ = WHSMask(freeze);
    std::cout << "ANNStrawHitUpdater " << " anncut " << mvacut_ << " null doca " << nulldoca_ << " freezeing " << freeze_ << std::endl;
    mva_->initMVA();
    mva_->showMVA();
    // outlier cuts: for now hardcode.  These should come from a central source FIXME
    mintdrift_ = -2.0;
    maxtdrift_ = 48.0;
    maxdoca_ = 5.0;
    maxresidpull_ = 10.0;
  }

  std::string const& ANNStrawHitUpdater::configDescription() {
    static std::string descrip( "Weight file, ANN cut, null hit doca, states to freeze");
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
      // outlier tests
      double rpull = (dinfo.driftDistance_ - fabs(tpdata.doca()))/sqrt(dinfo.driftDistanceVar() + tpdata.docaVar());
      if(tpdata.deltaT() < mintdrift_ ||
          tpdata.deltaT() > maxtdrift_ ||
          fabs(tpdata.doca()) > maxdoca_ ||
          abs(rpull) > maxresidpull_){
        whstate.state_ = WireHitState::inactive;
      } else {
        // invoke the ANN
        std::vector<Float_t> pars(9,0.0);
        // this order is given by the training
        pars[0] = fabs(tpdata.doca());
        pars[1] = dinfo.driftDistance_;
        pars[2] = chit.driftTime();
        pars[3] = 1000.0*chit.energyDep();
        pars[4] = rpull;
        pars[5] = sqrt(tpdata.docaVar());
        pars[6] = fabs(tpdata.dirDot());
        pars[7] = fabs(tpdata.particlePoca().Vect().Z());
        pars[8] = tpdata.particlePoca().Vect().Rho();
        float mvaout = mva_->evalMVA(pars);
        if(mvaout > mvacut_){
          whstate.state_ = tpdata.doca() > 0.0 ? WireHitState::right : WireHitState::left;
        } else {
          whstate.state_ = WireHitState::null;
        }
      }
      whstate.frozen_ = whstate.isIn(freeze_);
    }
    return whstate;
  }
}
