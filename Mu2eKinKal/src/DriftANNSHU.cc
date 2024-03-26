#include "Offline/Mu2eKinKal/inc/DriftANNSHU.hh"
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
    dtmvacut_ = std::get<4>(config);
    std::string freeze = std::get<5>(config);
    freeze_ = WHSMask(freeze);
    std::string flag = std::get<6>(config);
    flag_ = KKSHFlag(flag);
    diag_ = std::get<7>(config);
    if(diag_ > 0)
      std::cout << "DriftANNSHU LR sign weights " << std::get<0>(config) << " cut " << signmvacut_
        << " cluster weights " << std::get<0>(config) << " cut " << clustermvacut_ << " dt cut " << dtmvacut_
        << " freezing " << freeze_ << " flags " << flag << " diag " << diag_ << std::endl;
  }

  std::string const& DriftANNSHU::configDescription() {
    static std::string descrip( "LR ANN file, LR ANN cut, Drift ANN file, Drift ANN dresid cut, ANN dt cut, freeze states, flags, diag level");
    return descrip;
  }

  WireHitState DriftANNSHU::wireHitState(WireHitState const& input, ClosestApproachData const& tpdata, DriftInfo const& dinfo, ComboHit const& chit) const {
    WireHitState whstate = input;
    if(whstate.updateable(StrawHitUpdaters::DriftANN) && whstate.active()){
      // infer the ANN values
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
      if(diag_ > 2)std::cout << std::setw(8) << std::setprecision(5)
        << "Drift ANN inputs: doca, cdrift, sigdoca, TOTdrift, EDep "
          << spars[0] << " , "
          << spars[1] << " , "
          << spars[2] << " , "
          << spars[3] << " , "
          << spars[4] << " , "
          << " sign output " << signmvaout[0]
          << " drift output " << clustermvaout[0] << std::endl;
      whstate.quality_[WireHitState::sign] = signmvaout[0];
      whstate.quality_[WireHitState::drift] = clustermvaout[0];
      whstate.algo_  = StrawHitUpdaters::DriftANN;
      whstate.flag_ = flag_;
      bool setLR = clustermvaout[0] > clustermvacut_;
      bool annprob = flag_.hasAllProperties(KKSHFlag::annprob);
      if(!annprob){ // inteprept cut directly against the MVA output
        setLR &= signmvaout[0] > signmvacut_;
      } else {
        // interpret cut as scale for net benefit of LR assignment relative to null
        // Compute the expected variances for different scenarios
        double vr = dinfo.nullHitVar(); // variance assigned to a wire position constraint (null LR)
        double vs = dinfo.driftHitVar(); // variance used in weighting LR assigned hits
        double vx = tpdata.docavar_; // unsigned radius variance from other hits (existing fit)
        double vn = vx*vr/(vx+vr); // expected variance after assigning a null LR
        double vc = vx*vs/(vx+vs); // expected variance after assigning correct LR
        double ri = 2.0*dinfo.rDrift_*vx/(vs+vx); // weighted average unsigned radius from incorrect LR
        double vi = ri*ri + vc; // expected variance after assigning incorrect LR
        double relprob = 25*tan(1.5*signmvaout[0]); // this comes from a fit to p/(1-p): move this function to TrackerConditions TODO
        double lrprob = relprob/(1.0+relprob);// Probability LR assignment will be correct
        double vnet = lrprob*vc + (1.0-lrprob)*vi; // expected net radius variance if LR is assigned
        // use LR information if the (scaled) net variance assigning LR is smaller than the net variance adding a null hit
        setLR &= vnet*signmvacut_ < vn;
      }
      if(setLR){
        whstate.state_ = tpdata.doca() > 0.0 ? WireHitState::right : WireHitState::left;
        // only use dt constraint if the MVA are above the tighter cut.  This should be split for cluster, sign TODO
        if(clustermvaout[0] > dtmvacut_ && signmvaout[0] > dtmvacut_ ){
          whstate.flag_.merge(KKSHFlag::driftdt);
        } else {
          whstate.flag_.clear(KKSHFlag::driftdt);
        }
      } else {
        whstate.state_ = WireHitState::null;
      }
      whstate.frozen_ = whstate.isIn(freeze_);
      if (diag_ > 1)std::cout << "DriftANNSHU set hit " << whstate << std::endl;
    } else if (diag_ > 1) {
      std::cout << "DriftANNSHU skipping hit " << whstate << std::endl;
    }
    return whstate;
  }
}
