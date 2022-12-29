#include "Offline/Mu2eKinKal/inc/DriftANNSHU.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/Mu2eKinKal/inc/TrainDrift.hxx"
#include <cmath>
#include <array>

namespace mu2e {
  using KinKal::ClosestApproachData;
  using KinKal::VEC3;
  DriftANNSHU::DriftANNSHU(Config const& config) {
    ConfigFileLookupPolicy configFile;
    auto mvaWgtsFile = configFile(std::get<0>(config));
    mva_ = std::make_shared<TMVA_SOFIE_TrainDrift::Session>(mvaWgtsFile);
    mvacut_ = std::get<1>(config);
    std::string nulldvar = std::get<2>(config);
    nulldvar_ = WireHitState::nullDistVar(nulldvar);
    std::string allowed = std::get<3>(config);
    allowed_ = WHSMask(allowed);
    std::string freeze = std::get<4>(config);
    freeze_ = WHSMask(freeze);
    std::string totuse = std::get<5>(config);
    totuse_ = WireHitState::totUse(totuse);
    diag_ = std::get<6>(config);
    if(diag_ > 0)
      std::cout << "DriftANNSHU weights" << std::get<0>(config) << " cut " << mvacut_ << " null dist var " << nulldvar
        << " allowing " << allowed_ << " freezing " << freeze_ << " TOT use " << totuse << std::endl;
  }

  std::string const& DriftANNSHU::configDescription() {
    static std::string descrip( "Weight file, ANN cut, null hit doca, allowed states, TOT use, states to freeze, diag level");
    return descrip;
  }

  WireHitState DriftANNSHU::wireHitState(WireHitState const& input, ClosestApproachData const& tpdata, DriftInfo const& dinfo, ComboHit const& chit) const {
    WireHitState whstate = input;
    if(input.updateable(StrawHitUpdaters::DriftANN)){
      // invoke the ANN
      std::array<float,5> pars;
      // this order is given by the training
      pars[0] = fabs(tpdata.doca());
      pars[1] = dinfo.driftDistance_;
      pars[2] = tpdata.docaVar();
      pars[3] = chit.driftTime();
      pars[4] = chit.energyDep();
      auto mvaout = mva_->infer(pars.data());
      if(diag_ > 1){
        whstate.algo_  = StrawHitUpdaters::DriftANN;
        whstate.quality_ = mvaout[0];
      }
      if(mvaout[0] > mvacut_){
        if(allowed_.hasAnyProperty(WHSMask::drift)){
          whstate.state_ = tpdata.doca() > 0.0 ? WireHitState::right : WireHitState::left;
          whstate.algo_ = StrawHitUpdaters::DriftANN;
          whstate.totuse_ = totuse_;
          whstate.quality_ = mvaout[0];
       }
      } else {
        if(allowed_.hasAnyProperty(WHSMask::null)){
          whstate.state_ = WireHitState::null;
          whstate.algo_ = StrawHitUpdaters::DriftANN;
          whstate.totuse_ = totuse_;
          whstate.nulldvar_ = nulldvar_;
          whstate.quality_ = mvaout[0];
        }
      }
      if(whstate.algo_ == StrawHitUpdaters::DriftANN)whstate.frozen_ = whstate.isIn(freeze_);
    }
    return whstate;
  }
}
