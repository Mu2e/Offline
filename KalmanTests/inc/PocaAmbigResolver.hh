//
// resolver that sets every hit ambiguity according to POCA on each iteration
// $Id: PocaAmbigResolver.hh,v 1.1 2012/05/14 19:20:02 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/05/14 19:20:02 $
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
      virtual void resolveTrk(TrkKalFit& kfit) const;
    private:
  };
}
#endif
