//
// base class to resolve hit ambiguities.  This class also interacts with the
// simulated annealing.
//
// Original author: David Brown (LBNL), 2012
//
// $Id: AmbigResolver.hh,v 1.4 2014/08/01 18:56:10 gandr Exp $
// $Author: gandr $ 
// $Date: 2014/08/01 18:56:10 $
//
#ifndef AmbigResolver_HH
#define AmbigResolver_HH
#include "BTrk/BaBar/BaBar.hh"

class KalRep;
class TrkSimpTraj;
namespace mu2e {
  class TrkStrawHit;

  class AmbigResolver {
    public:
    explicit AmbigResolver(double extErr);
    virtual ~AmbigResolver() = 0;

// resolve a track.  Depending on the configuration, this might
// update the hit state and the t0 value.

    virtual void resolveTrk(KalRep* kfit) const = 0;

    protected:
// init hit external errors for simulated annealing
    virtual void initHitErrors(KalRep* kfit) const ;
// find the local trajectory piece computed from the fit excluding a particular set of hits.
// the hits are assumed to be contiguous
    const TrkSimpTraj* findTraj(std::vector<TrkStrawHit*> const& phits, const KalRep* krep) const;
    double _extErr; // external hit error
  };
}

#endif
