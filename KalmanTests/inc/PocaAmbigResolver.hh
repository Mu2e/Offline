//
// resolver that sets every hit ambiguity according to POCA on each iteration
// $Id: PocaAmbigResolver.hh,v 1.2 2012/08/31 22:39:00 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/08/31 22:39:00 $
//
#ifndef PocaAmbigResolver_HH
#define PocaAmbigResolver_HH
#include "BaBar/BaBar.hh"
#include "KalmanTests/inc/AmbigResolver.hh"

namespace mu2e {

  class PocaAmbigResolver : public AmbigResolver {
    public:
// construct from parameter set
      explicit PocaAmbigResolver(fhicl::ParameterSet const& pset);
      virtual ~PocaAmbigResolver();
      virtual void resolveTrk(KalFitResult& kfit) const;
    private:
  };
}
#endif
