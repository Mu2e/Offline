#include "Offline/Mu2eKinKal/inc/CombinatoricStrawHitUpdater.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include <algorithm>
#include <iostream>
namespace mu2e {
  NullHitInfo CombinatoricStrawHitUpdater::nullHitInfo(StrawResponse const& sresponse, Straw const& straw) const {
    NullHitInfo nhinfo;
    // compute time and distance parameters used for null ambiguity (wire constraint)
    double vdriftinst = sresponse.driftInstantSpeed(straw.id(),nulldoca_,0.0,true);
    nhinfo.toff_ = 0.5*nulldoca_/vdriftinst;
    nhinfo.tvar_ = 0.25*dvar_/(vdriftinst*vdriftinst);
    nhinfo.dvar_ = dvar_;
    nhinfo.usetime_ = nulltime_;
    return nhinfo;
  }


  double CombinatoricStrawHitUpdater::penalty(WireHitState const& whstate) const {
    if(!whstate.active())
      return inactivep_;
    else if(!whstate.useDrift())
      return nullp_;
    else
      return 0.0;
  }

  ClusterScore CombinatoricStrawHitUpdater::selectBest(ClusterScoreCOL& cscores) const {
    // sort the results by chisquared
    std::sort(cscores.begin(),cscores.end(), ClusterScoreComp());
    auto best = cscores.front();
    // set the algorithm
    for(auto& whs : best.hitstates_)whs.algo_ = StrawHitUpdaters::Combinatoric;
    // merge the configurations that have nearly degenerate scores; this takes the most conservative option
    auto test=cscores.begin(); ++test;
    while(test != cscores.end() && test->chi2_.chisqPerNDOF() - cscores.front().chi2_.chisqPerNDOF() < minDeltaChi2()){
      best.merge(*test);
      ++test;
    }
    // optionally freeze unambiguous states
    if(freeze_){
      // look for pairs of unambiguous (opposite drift) hits, and freeze their state
      for( size_t ihit=0;ihit < best.hitstates_.size();++ihit){
        auto& hit1 = best.hitstates_[ihit];
        for( size_t jhit=ihit+1;jhit < best.hitstates_.size();++jhit){
          auto& hit2 = best.hitstates_[jhit];
          if(hit1.useDrift() && hit2.useDrift() && hit1 != hit2){
            hit1.frozen_ = true;
            hit2.frozen_ = true;
          }
        }
      }
      //      for(auto& whs : best.hitstates_) {
      //        whs.frozen_ = whs.useDrift();
      //      }
    }
    return best;
  }

  void ClusterScore::merge(ClusterScore const& other) {
    for( size_t ihit=0;ihit < hitstates_.size();++ihit){
      auto& myhit = hitstates_[ihit];
      auto const& otherhit = other.hitstates_[ihit];
      if(myhit != otherhit){
        if(myhit.isInactive() || otherhit.isInactive()) // one is inactive: deactivate
          myhit.state_ = WireHitState::inactive;
        else // every other case, set to null
          myhit.state_ = WireHitState::null;
      }
    }
  }

  std::ostream& operator <<(std::ostream& os, ClusterScore const& cscore ) {
    os << "ClusterScore " << cscore.chi2_ << " states: ";
    for(auto whstate : cscore.hitstates_) std::cout << "  " << whstate.state_;
    std::cout << std::endl;
    return os;
  }

}
