//
// dummy 'resolver' that leaves every hit ambiguity unchanged
//
// Original author: David Brown (LBNL), 2012
//
// $Id: FixedAmbigResolver.hh,v 1.3 2014/08/01 18:56:10 gandr Exp $
// $Author: gandr $ 
// $Date: 2014/08/01 18:56:10 $
//
#ifndef FixedAmbigResolver_HH
#define FixedAmbigResolver_HH
#include "BTrk/BaBar/BaBar.hh"
#include "TrkReco/inc/AmbigResolver.hh"
#ifndef __GCCXML__
#include "fhiclcpp/ParameterSet.h"
#endif/*__GCCXML__*/

namespace mu2e {

  class FixedAmbigResolver : public AmbigResolver {
    public:
// construct from parameter set
#ifndef __GCCXML__
    explicit FixedAmbigResolver(fhicl::ParameterSet const& pset,double tmpErr);
#endif/*__GCCXML__*/
      virtual ~FixedAmbigResolver();
    virtual bool resolveTrk(KalRep* kfit) const;
    private:
      bool _neutralize; // if true, set the initial ambiguity to 0
  };
}
#endif
