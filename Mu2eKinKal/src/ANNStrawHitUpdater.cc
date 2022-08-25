#include "Offline/Mu2eKinKal/inc/ANNStrawHitUpdater.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include <cmath>

namespace mu2e {
  using KinKal::ClosestApproachData;
  using KinKal::VEC3;
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
      if(mvaout > mvacut_){
        whstate.state_ = tpdata.doca() > 0.0 ? WireHitState::right : WireHitState::left;
        whstate.frozen_ = freeze_;
      } else {
        whstate.state_ = WireHitState::null;
        whstate.nhmode_ = nhmode_;
        // compute time and distance parameters used for null ambiguity (wire constraint)
        if(nhmode_ == WireHitState::combo){
          whstate.dvar_ =  2.0833; // (2*rstraw)^2/12   Should come from TrackerGeom TODO
        } else {
          whstate.dvar_ = dvar_;
        }
      }
    }
    return whstate;
  }
}
