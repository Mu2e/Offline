//
// class to resolve hit ambiguities by a single hit, assuming a reasonable track
// fit as input
// 
// Original author: David Brown (LBNL), 2012
//
// $Id: HitAmbigResolver.hh,v 1.4 2014/08/01 18:56:10 gandr Exp $
// $Author: gandr $ 
// $Date: 2014/08/01 18:56:10 $
//
#ifndef HitAmbigResolver_HH
#define HitAmbigResolver_HH
#include "BTrk/BaBar/BaBar.hh"
#include "TrkReco/inc/AmbigResolver.hh"
#ifndef __GCCXML__
#include "fhiclcpp/ParameterSet.h"
#endif/*__GCCXML__*/
#include <cstddef>
#include <vector>

class TrkDifTraj;
class KalRep;

namespace mu2e {

  class HitAmbigResolver : public AmbigResolver {
    public:
      enum trajtype {reftraj=0};
// construct from parameter set
#ifndef __GCCXML__
    explicit HitAmbigResolver(fhicl::ParameterSet const&, double tmpErr);
#endif/*__GCCXML__*/
      virtual ~HitAmbigResolver();
// resolve a track.  Depending on the configuration, this might
// update the hit state and the t0 value.
    virtual bool resolveTrk(KalRep* kfit) const;
    private:
// penalty function depends on the drift radius
      double penaltyError(double rdrift) const;
      double _mindrift; // minimum drift to assign an ambiguity.  Below this, an ambiguity of '0' is defined
      double _zeropenalty; // special penalty for drifts below the minimum
      bool _penalty; // apply penalty or notA
// exponential + linear fit to ambiguity mis-assignment
      double _expnorm;
      double _lambda;
      double _offset;
      double _slope;
      int _debug;
  };
}
#endif
