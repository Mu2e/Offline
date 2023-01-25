#include "Offline/Mu2eKinKal/inc/Chi2SHU.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include <algorithm>
#include <iostream>
namespace mu2e {

  Chi2SHU::Chi2SHU(Config const& config) {
    csize_ = std::get<0>(config);
    inactivep_ = std::get<1>(config);
    nullp_ = std::get<2>(config);
    mindchi2_ = std::get<3>(config);
    std::string nulldvar = std::get<4>(config);
    nulldvar_ = WireHitState::nullDistVar(nulldvar);
    std::string states = std::get<5>(config);
    WHSMask allowed(states);
    std::string freeze = std::get<6>(config);
    freeze_ = WHSMask(freeze);
    std::string unfreeze = std::get<7>(config);
    unfreeze_ = WHSMask(unfreeze);
    diag_ = std::get<8>(config);
    if(allowed.hasAnyProperty(WHSMask::inactive)) allowed_.emplace_back(WireHitState::inactive,StrawHitUpdaters::Chi2,nulldvar_);
    if(allowed.hasAnyProperty(WHSMask::null)) allowed_.emplace_back(WireHitState::null,StrawHitUpdaters::Chi2,nulldvar_);
    if(allowed.hasAnyProperty(WHSMask::drift)){
      allowed_.emplace_back(WireHitState::left,StrawHitUpdaters::Chi2,nulldvar_);
      allowed_.emplace_back(WireHitState::right,StrawHitUpdaters::Chi2,nulldvar_);
    }
    if(diag_ > 0)std::cout << "Chi2SHU, inactive penalty " << inactivep_ << " null penalty " << nullp_ << " min dchi2 " << mindchi2_ << " null dist var " << nulldvar << " allowed states" << allowed << " states to freeze " << freeze_  << " states to unfreeze" << unfreeze_ << std::endl;
  }

  // set the state of unambiguous hits to their drift value.
  ClusterState Chi2SHU::selectBest(ClusterStateCOL& cstates) const {
    // sort the results by chisquared
    std::sort(cstates.begin(),cstates.end(), ClusterStateComp());

    auto best = cstates.front();
    // define 'quality' as the chi2 gap between the 2 best states
    double quality(-1.0);
    if(cstates.size() > 1)quality = cstates[1].chi2_.chisq() - cstates[0].chi2_.chisq();

    // merge the configurations that have nearly degenerate states; this takes the most conservative option
    auto test=cstates.begin(); ++test;
    while(test != cstates.end() && test->chi2_.chisq() - cstates.front().chi2_.chisq() < mindchi2_){
      best.merge(*test);
      ++test;
    }
    for(auto& whs : best.hitstates_){
      whs.frozen_ =  whs.isIn(freeze_);
      whs.quality_[WireHitState::chi2] = quality;
    }
    if(diag_ > 1){
      std::cout << "Best Cluster " << best.chi2_ << " hit states ";
      for(auto whs : best.hitstates_)std::cout << "  " << whs.state_;
      std::cout << std::endl;
    }
    if(diag_ > 2){
      for(auto const&  cstate: cstates) {
        std::cout << "Combi " << cstate.chi2_ << " hit states ";
        for(auto whstate : cstate.hitstates_) std::cout << "  " << whstate.state_;
        std::cout << std::endl;
      }
    }
    // define the quality as the difference between the best and the next-best assignment
    return best;
  }

  void ClusterState::merge(ClusterState const& other) {
    for( size_t ihit=0;ihit < hitstates_.size();++ihit){
      auto& mystate = hitstates_[ihit];
      auto const& otherstate = other.hitstates_[ihit];
      if(mystate != otherstate){
        if(mystate.isInactive() || otherstate.isInactive()) // one is inactive: deactivate
          mystate.state_ = WireHitState::inactive;
        else // every other case, set to null
          mystate.state_ = WireHitState::null;
      }
    }
  }

  std::ostream& operator <<(std::ostream& os, ClusterState const& cstate ) {
    os << "ClusterState " << cstate.chi2_ << " states: ";
    for(auto whstate : cstate.hitstates_) std::cout << "  " << whstate.state_;
    std::cout << std::endl;
    return os;
  }

  std::string const& Chi2SHU::configDescription() {
    static std::string descrip("Min Cluster Size, Inactive hit x^2 penalty, Null ambiguity x^2 penalty, Minimum significant x^2 difference, minimum drift DOCA, allowed states, states to freeze, states to unfreeze, diag level");

    return descrip;
  }

}
