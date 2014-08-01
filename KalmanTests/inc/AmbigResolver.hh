//
// base class to resolve hit ambiguities 
//
// $Id: AmbigResolver.hh,v 1.4 2014/08/01 18:56:10 gandr Exp $
// $Author: gandr $ 
// $Date: 2014/08/01 18:56:10 $
//
#ifndef AmbigResolver_HH
#define AmbigResolver_HH
#include "BaBar/BaBar.hh"
#ifndef __GCCXML__
#include "fhiclcpp/ParameterSet.h"
#endif/*__GCCXML__*/

class KalRep;
class TrkSimpTraj;
namespace mu2e {
  class KalFitResult;
  class TrkStrawHit;

  class AmbigResolver {
    public:
// construct from parameter set
#ifndef __GCCXML__
      explicit AmbigResolver(fhicl::ParameterSet const&);
#endif/*__GCCXML__*/
      virtual ~AmbigResolver() = 0;
// resolve a track.  Depending on the configuration, this might
// update the hit state and the t0 value.
      virtual void resolveTrk(KalFitResult& kfit) const = 0;
    protected:
// find the local trajectory piece computed from the fit excluding a particular set of hits.
// the hits are assumed to be contiguous
      const TrkSimpTraj* findTraj(std::vector<TrkStrawHit*> const& phits, const KalRep* krep) const;
  };
}

#endif
