#include "Offline/Mu2eKinKal/inc/ANNStrawHitUpdater.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include <cmath>

namespace mu2e {
  using KinKal::ClosestApproachData;
  using KinKal::VEC3;
  WireHitState ANNStrawHitUpdater::wireHitState(ClosestApproachData const& tpdata, DriftInfo const& dinfo, ComboHit const& chit) const {
    std::vector<Float_t> pars(8,0.0);
    // this order is given by the training
    pars[0] = abs(tpdata.doca());
    pars[1] = dinfo.driftDistance_;
    pars[2] = chit.driftTime();
    pars[3] = 1000.0*chit.energyDep();
    pars[4] = tpdata.docaVar();
    pars[5] = tpdata.dirDot();
    pars[6] = tpdata.particlePoca().Vect().Rho();
    pars[7] = tpdata.particlePoca().Vect().Z();
    float mvaout = mva_->evalMVA(pars);
    WireHitState whstate(WireHitState::inactive,StrawHitUpdaters::ANN);
    if(mvaout > mvacut_){
      whstate.state_ = tpdata.doca() > 0.0 ? WireHitState::right : WireHitState::left;
    } else { //  Need a way to decide to disable the hit as well TODO
      whstate.state_ = WireHitState::null;
    }
    // compute time and distance parameters used for null ambiguity (wire constraint)
    auto& nhinfo = whstate.nhinfo_;
    nhinfo.tmode_ = nhtmode_;
    if(nhtmode_ == NullHitInfo::usecombo){
      nhinfo.toff_ =  0.0;
      nhinfo.tvar_ = 50.0; // should come from ComboHit FIXME
      nhinfo.dvar_ = dvar_;
    } else if(nhtmode_ == NullHitInfo::usedoca){
      double vdrift = dinfo.driftVelocity_;
      nhinfo.toff_ = std::max(0.0,dinfo.driftDistance_)/vdrift; // this calculation is unreliable currently
      nhinfo.dvar_ = dvar_;
      nhinfo.tvar_ = dvar_/(vdrift*vdrift);
    }
    return whstate;
  }
}
