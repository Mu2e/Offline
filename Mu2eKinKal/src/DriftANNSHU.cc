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
    clusterdtmvacut_ = std::get<4>(config);
    std::string nulldvar = std::get<5>(config);
    nulldvar_ = WireHitState::nullDistVar(nulldvar);
    std::string allowed = std::get<6>(config);
    allowed_ = WHSMask(allowed);
    std::string freeze = std::get<7>(config);
    freeze_ = WHSMask(freeze);
    std::string totuse = std::get<8>(config);
    totuse_ = WireHitState::totUse(totuse);
    diag_ = std::get<9>(config);
    if(diag_ > 0)
      std::cout << "DriftANNSHU LR sign weights " << std::get<0>(config) << " cut " << signmvacut_
        << " cluster weights " << std::get<0>(config) << " cuts " << clustermvacut_ << " , " << clusterdtmvacut_
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
      spars[1] = dinfo.cDrift_;
      spars[2] = sqrt(std::max(0.0,tpdata.docaVar()));
      spars[3] = chit.driftTime();
      // normalize edep.
      // For sign, noralize only to the crossing angle, as there it serves as an estimate of the drift radius
      double sint = sqrt(1.0-tpdata.dirDot()*tpdata.dirDot());
      spars[4] = chit.energyDep()*sint;
      auto signmvaout = signmva_->infer(spars.data());
      cpars[0] = fabs(tpdata.doca());
      cpars[1] = dinfo.cDrift_;
      cpars[2] = chit.driftTime();
      // For drift quality, normalize to the estimated path length through the straw, as that measures the clustering effects
      double plen = sqrt(std::max(0.25, 6.25-dinfo.rDrift_*dinfo.rDrift_))/sint;
      cpars[3] = chit.energyDep()/plen;
      auto clustermvaout = clustermva_->infer(cpars.data());
      whstate.quality_[WireHitState::sign] = signmvaout[0];
      whstate.quality_[WireHitState::drift] = clustermvaout[0];
      if(diag_ > 1) whstate.algo_  = StrawHitUpdaters::DriftANN;
      // Compute the cost to benefit ratio of correct vs incorrect LR assignment, and only sign if the propbability of a wrong
      // assignment (estimated by the MVA) is less than that (times a scale factor)
      double vr = dinfo.nullHitVar();
      double vh = dinfo.unsignedDriftVar();
      double vx = tpdata.docavar_;
      double vn = vx*vr/(vx+vr); // expected variance assigning a null LR
      double vc = vx*vh/(vx+vh); // expected variance assigning correct LR
      double vi = 4.0*vr*vc/(vx+vh); // expected variance assigning incorrect LR
      double signcost = (vi - vn)/(vn - vc); // net
//      double signcost = dinfo.nullHitVar()/tpdata.docavar_;
      // the following should be a calibration function TODO
      double signprob = 25*tan(1.5*signmvaout[0]); // this is a fit to p/(1-p)
//      if(signmvaout[0] > signmvacut_ && clustermvaout[0] > clustermvacut_){
      if(signprob > signcost*signmvacut_ && clustermvaout[0] > clustermvacut_){
        if(allowed_.hasAnyProperty(WHSMask::drift)){
          if(clustermvaout[0] > clusterdtmvacut_)
            whstate.state_ = tpdata.doca() > 0.0 ? WireHitState::rightdt : WireHitState::leftdt;
          else
            whstate.state_ = tpdata.doca() > 0.0 ? WireHitState::right : WireHitState::left;
          whstate.algo_ = StrawHitUpdaters::DriftANN;
          whstate.totuse_ = totuse_;
          whstate.nulldvar_ = nulldvar_;
        }
      } else {
        if(allowed_.hasAnyProperty(WHSMask::null)){
          whstate.state_ = WireHitState::null;
          whstate.algo_ = StrawHitUpdaters::DriftANN;
          whstate.totuse_ = totuse_;
          whstate.nulldvar_ = nulldvar_;
        }
      }
      if(whstate.algo_ == StrawHitUpdaters::DriftANN)whstate.frozen_ = whstate.isIn(freeze_);
    }
    return whstate;
  }
}
