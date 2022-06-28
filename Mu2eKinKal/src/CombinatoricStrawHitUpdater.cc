#include "Offline/Mu2eKinKal/inc/CombinatoricStrawHitUpdater.hh"
#include <algorithm>
namespace mu2e {

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
    // Test if the best solutions are nearly degenerate, choose the most conservative option
    auto best = cscores.begin();
    double bestrank = wireHitRank(best->hitstates_);

    auto test=best; ++test;
    while(test->chi2_.chisqPerNDOF() - cscores.front().chi2_.chisqPerNDOF() < minDeltaChi2() && test != cscores.end()){
      double testrank = wireHitRank(test->hitstates_);
      if(testrank > bestrank){
        best = test;
        bestrank = testrank;
      }
      ++test;
    }
    return *best;
 }

  double CombinatoricStrawHitUpdater::wireHitRank(WHSCOL const& hitstates) const {
    double rank(0.0);
    for(auto const& whs : hitstates){
      switch (whs.state_) {
        case WireHitState::inactive:
          rank += 2.0;
          break;
        case WireHitState::null:
          rank += 1.0;
          break;
        default:
          break;
      }
    }
    return rank;
  }
}
