//
// base class to resolve hit ambiguities.  This class also interacts with the
// simulated annealing.
//
// Original author: David Brown (LBNL), 2012
//
//
#ifndef AmbigResolver_HH
#define AmbigResolver_HH
// #include "BTrk/BaBar/BaBar.hh"

#include <vector>

class KalRep;
class TrkSimpTraj;
namespace mu2e {
  class TrkStrawHit;
  class AmbigResolver {
    public:
      explicit AmbigResolver(double tmpErr);
      virtual ~AmbigResolver() = 0;

      // resolve a track.  Depending on the configuration, this might
      // update the hit state and the t0 value.

      virtual bool resolveTrk(KalRep* kfit) const = 0;

    protected:
      // reset penalty errors
      virtual void initHitErrors(KalRep* kfit) const ;
      // find the local trajectory piece computed from the fit excluding a particular set of hits.
      // the hits are assumed to be contiguous
      const TrkSimpTraj* findTraj(std::vector<TrkStrawHit*> const& phits, const KalRep* krep) const;
      double _tmpErr; // hit error associated with annealing 'temperature'
  };
}

#endif
