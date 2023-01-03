#include "Offline/Mu2eKinKal/inc/DriftANNSHU.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include <cmath>
#include <array>

namespace mu2e {
  using KinKal::ClosestApproachData;
  using KinKal::VEC3;
  DriftANNSHU::DriftANNSHU(Config const& config) {
    ConfigFileLookupPolicy configFile;
    auto signmvaWgtsFile = configFile(std::get<0>(config));
    signmva_ = std::make_shared<TMVA_SOFIE_TrainSign::Session>(signmvaWgtsFile);
    signmvacut_ = std::get<1>(config);
    auto clustermvaWgtsFile = configFile(std::get<2>(config));
    clustermva_ = std::make_shared<TMVA_SOFIE_TrainCluster::Session>(clustermvaWgtsFile);
    clustermvacut_ = std::get<3>(config);
    std::string nulldvar = std::get<4>(config);
    nulldvar_ = WireHitState::nullDistVar(nulldvar);
    std::string allowed = std::get<5>(config);
    allowed_ = WHSMask(allowed);
    std::string freeze = std::get<6>(config);
    freeze_ = WHSMask(freeze);
    std::string totuse = std::get<7>(config);
    totuse_ = WireHitState::totUse(totuse);
    diag_ = std::get<8>(config);
    if(diag_ > 0)
      std::cout << "DriftANNSHU LR sign weights " << std::get<0>(config) << " cut " << signmvacut_
        << " cluster weights " << std::get<0>(config) << " cut " << clustermvacut_
        << " null dist var " << nulldvar
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
      std::array<float,5> spars;
      std::array<float,4> cpars;
      // this order is given by the training
      spars[0] = fabs(tpdata.doca());
      spars[1] = dinfo.driftDistance_;
      spars[2] = sqrt(std::max(0.0,tpdata.docaVar()));
      spars[3] = chit.driftTime();
      spars[4] = chit.energyDep();
      auto signmvaout = signmva_->infer(spars.data());
      cpars[0] = fabs(tpdata.doca());
      cpars[1] = dinfo.driftDistance_;
      cpars[2] = chit.driftTime();
      cpars[3] = chit.energyDep();
     auto clustermvaout = clustermva_->infer(cpars.data());
      if(diag_ > 1){
        whstate.algo_  = StrawHitUpdaters::DriftANN;
        whstate.quality_ = signmvaout[0]; // need an array here TODO
      }
      if(signmvaout[0] > signmvacut_ && clustermvaout[0] > clustermvacut_){
        if(allowed_.hasAnyProperty(WHSMask::drift)){
          whstate.state_ = tpdata.doca() > 0.0 ? WireHitState::right : WireHitState::left;
          whstate.algo_ = StrawHitUpdaters::DriftANN;
          whstate.totuse_ = totuse_;
          whstate.quality_ = signmvaout[0];
       }
      } else {
        if(allowed_.hasAnyProperty(WHSMask::null)){
          whstate.state_ = WireHitState::null;
          whstate.algo_ = StrawHitUpdaters::DriftANN;
          whstate.totuse_ = totuse_;
          whstate.nulldvar_ = nulldvar_;
          whstate.quality_ = signmvaout[0];
        }
      }
      if(whstate.algo_ == StrawHitUpdaters::DriftANN)whstate.frozen_ = whstate.isIn(freeze_);
    }
    return whstate;
  }
}
