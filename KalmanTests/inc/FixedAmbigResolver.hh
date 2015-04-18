//
// dummy 'resolver' that leaves every hit ambiguity unchanged
// $Id: FixedAmbigResolver.hh,v 1.3 2014/08/01 18:56:10 gandr Exp $
// $Author: gandr $ 
// $Date: 2014/08/01 18:56:10 $
//
#ifndef FixedAmbigResolver_HH
#define FixedAmbigResolver_HH
#include "BaBar/BaBar.hh"
#include "KalmanTests/inc/AmbigResolver.hh"

namespace mu2e {

  class FixedAmbigResolver : public AmbigResolver {
    public:
// construct from parameter set
#ifndef __GCCXML__
    explicit FixedAmbigResolver(fhicl::ParameterSet const& pset, double ExtErr, int Iter);
#endif/*__GCCXML__*/
      virtual ~FixedAmbigResolver();
    virtual void resolveTrk(KalFitResult& kfit, int Final) const;
    private:
      bool _neutralize; // if true, set the initial ambiguity to 0
  };
}
#endif
