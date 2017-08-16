//
// resolver that sets every hit ambiguity according to POCA on each iteration
//
// Original author: David Brown (LBNL), 2012
//
// $Id: PocaAmbigResolver.hh,v 1.3 2014/08/01 18:56:10 gandr Exp $
// $Author: gandr $ 
// $Date: 2014/08/01 18:56:10 $
//
#ifndef PocaAmbigResolver_HH
#define PocaAmbigResolver_HH
#include "BTrk/BaBar/BaBar.hh"
#include "TrkReco/inc/AmbigResolver.hh"
#ifndef __GCCXML__
#include "fhiclcpp/ParameterSet.h"
#endif/*__GCCXML__*/

namespace mu2e {

  class PocaAmbigResolver : public AmbigResolver {
  public:
    // construct from parameter set
#ifndef __GCCXML__
    explicit PocaAmbigResolver(fhicl::ParameterSet const& pset);
#endif/*__GCCXML__*/
    virtual ~PocaAmbigResolver();
    virtual bool resolveTrk(KalRep* kfit) const;
  private:
  };
}
#endif
