//
// base class to resolve hit ambiguities 
//
// $Id: AmbigResolver.hh,v 1.2 2012/07/25 20:56:57 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/07/25 20:56:57 $
//
#ifndef AmbigResolver_HH
#define AmbigResolver_HH
#include "BaBar/BaBar.hh"
#include "fhiclcpp/ParameterSet.h"

class KalRep;
class TrkSimpTraj;
namespace mu2e {
  class TrkKalFit;
  class TrkStrawHit;

  class AmbigResolver {
    public:
// construct from parameter set
      explicit AmbigResolver(fhicl::ParameterSet const&);
      virtual ~AmbigResolver() = 0;
// resolve a track.  Depending on the configuration, this might
// update the hit state and the t0 value.
      virtual void resolveTrk(TrkKalFit& kfit) const = 0;
    protected:
// find the local trajectory piece computed from the fit excluding a particular set of hits.
// the hits are assumed to be contiguous
      const TrkSimpTraj* findTraj(std::vector<TrkStrawHit*> const& phits, const KalRep* krep) const;
  };
}

#endif
