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
    // pick the best configuration and update the hits
    // TODO: Test if the best solutions are nearly degenerate, and if so and they involve flipping a single hit, chose the most conservative option
    // This is especially important for single-hit 'clusters'
    return cscores.front();
  }
}
