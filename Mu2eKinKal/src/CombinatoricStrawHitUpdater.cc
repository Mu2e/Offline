#include "Offline/Mu2eKinKal/inc/CombinatoricStrawHitUpdater.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include <algorithm>
#include <iostream>
namespace mu2e {

  CombinatoricStrawHitUpdater::CombinatoricStrawHitUpdater(CSHUConfig const& cshuconfig) {
    csize_ = std::get<0>(cshuconfig);
    inactivep_ = std::get<1>(cshuconfig);
    nullp_ = std::get<2>(cshuconfig);
    mindchi2_ = std::get<3>(cshuconfig);
    nulldoca_ = std::get<4>(cshuconfig);
    std::string states = std::get<5>(cshuconfig);
    auto allowed = WHSMask(states);
    std::string freeze = std::get<6>(cshuconfig);
    freeze_ = WHSMask(freeze);
    unfreeze_ = std::get<7>(cshuconfig);
    diag_ = std::get<8>(cshuconfig);
    if(allowed.hasAnyProperty(WHSMask::inactive)) allowed_.emplace_back(WireHitState::inactive,StrawHitUpdaters::Combinatoric,nulldoca_);
    if(allowed.hasAnyProperty(WHSMask::null)) allowed_.emplace_back(WireHitState::null,StrawHitUpdaters::Combinatoric,nulldoca_);
    if(allowed.hasAnyProperty(WHSMask::drift)){
      allowed_.emplace_back(WireHitState::left,StrawHitUpdaters::Combinatoric,nulldoca_);
      allowed_.emplace_back(WireHitState::right,StrawHitUpdaters::Combinatoric,nulldoca_);
    }
    if(diag_ > 0)std::cout << "CombinatoricStrawHitUpdater " << inactivep_ << " " << nullp_ << " " << mindchi2_ << " " << nulldoca_ << " allowed states" << allowed << " states to freeze " << freeze_  << " unfreeze input? " << unfreeze_ << std::endl;
  }

  // set the state of unambiguous hits to their drift value.
  ClusterState CombinatoricStrawHitUpdater::selectBest(ClusterStateCOL& cstates) const {
    // sort the results by chisquared
    std::sort(cstates.begin(),cstates.end(), ClusterStateComp());

    auto best = cstates.front();
    double quality(-1.0);
    if(cstates.size() > 1)quality = cstates[1].chi2_.chisq() - cstates[0].chi2_.chisq();

    // merge the configurations that have nearly degenerate states; this takes the most conservative option
    auto test=cstates.begin(); ++test;
    while(test != cstates.end() && test->chi2_.chisq() - cstates.front().chi2_.chisq() < minDeltaChi2()){
      best.merge(*test);
      ++test;
    }
    for(auto& whs : best.hitstates_){
      whs.frozen_ =  whs.isIn(freeze_);
      whs.quality_ = quality;
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

  std::string const& CombinatoricStrawHitUpdater::configDescription() {
    static std::string descrip("Min Cluster Size, Inactive hit x^2 penalty, Null ambiguity x^2 penalty, Minimum significant x^2 difference, minimum drift DOCA, allowed states, states to freeze, unfreeze input, diag level");

    return descrip;
  }

}
