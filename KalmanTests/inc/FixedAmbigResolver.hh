//
// dummy 'resolver' that leaves every hit ambiguity unchanged
// $Id: FixedAmbigResolver.hh,v 1.2 2012/08/31 22:38:59 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/08/31 22:38:59 $
//
#ifndef FixedAmbigResolver_HH
#define FixedAmbigResolver_HH
#include "BaBar/BaBar.hh"
#include "KalmanTests/inc/AmbigResolver.hh"

namespace mu2e {

  class FixedAmbigResolver : public AmbigResolver {
    public:
// construct from parameter set
      explicit FixedAmbigResolver(fhicl::ParameterSet const& pset);
      virtual ~FixedAmbigResolver();
      virtual void resolveTrk(KalFitResult& kfit) const;
    private:
      bool _neutralize; // if true, set the initial ambiguity to 0
  };
}
#endif
