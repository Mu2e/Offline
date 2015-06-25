//
// resolver that sets every hit ambiguity according to POCA on each iteration
// $Id: PocaAmbigResolver.hh,v 1.3 2014/08/01 18:56:10 gandr Exp $
// $Author: gandr $ 
// $Date: 2014/08/01 18:56:10 $
//
#ifndef PocaAmbigResolver_HH
#define PocaAmbigResolver_HH
#include "BTrk/BaBar/BaBar.hh"
#include "KalmanTests/inc/AmbigResolver.hh"

namespace mu2e {

  class PocaAmbigResolver : public AmbigResolver {
  public:
    // construct from parameter set
#ifndef __GCCXML__
    explicit PocaAmbigResolver(fhicl::ParameterSet const& pset, double ExtErr, int Iter);
#endif/*__GCCXML__*/
    virtual ~PocaAmbigResolver();
    virtual void resolveTrk(KalFitResult& kfit, int Final) const;
  private:
  };
}
#endif
