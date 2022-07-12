#include "Offline/Mu2eKinKal/inc/CombinatoricStrawHitUpdater.hh"
#include <algorithm>
#include <iostream>
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
    auto best = cscores.front();
    double bestrank = wireHitRank(best.hitstates_);

    auto test=cscores.begin(); ++test;
    while(test != cscores.end() &&
        test->chi2_.chisqPerNDOF() - cscores.front().chi2_.chisqPerNDOF() < minDeltaChi2()){
      if(allownull_){
        double testrank = wireHitRank(test->hitstates_);
        if(testrank > bestrank){
          if(diag_ > 0) std::cout << "Replacing " << best << " with " << *test << std::endl;
          best = *test;
          bestrank = testrank;
        }
        ++test;
      } else {
        // if we didn't explicitly assign null hits, assign them now if the difference between left and right is small
        while(test != cscores.end() &&
            test->chi2_.chisqPerNDOF() - cscores.front().chi2_.chisqPerNDOF() < minDeltaChi2()){
          for(size_t ihit=0;ihit<best.hitstates_.size();++ihit){
            // if this hit flips without appreciably changing chisquared, set it null
            if(best.hitstates_[ihit].state_*test->hitstates_[ihit].state_ == WireHitState::left){
              if(diag_ > 0) std::cout << "Nulling hit  " << ihit
                << " dchisq " << test->chi2_.chisqPerNDOF() - cscores.front().chi2_.chisqPerNDOF() <<std::endl;
              best.hitstates_[ihit] = WireHitState::null;
            }
          }
        }
      }
    }
    return best;
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

  std::ostream& operator <<(std::ostream& os, ClusterScore const& cscore ) {
    os << "ClusterScore " << cscore.chi2_ << " states: ";
    for(auto whstate : cscore.hitstates_) std::cout << "  " << whstate.state_;
    std::cout << std::endl;
    return os;
  }


}
