//
// class to resolve hit ambiguities one hit at a time, assuming a reasonable track
// fit as input
//
// $Id: FixedAmbigResolver.cc,v 1.3 2012/08/31 22:39:00 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/08/31 22:39:00 $
//
#include "KalmanTests/inc/FixedAmbigResolver.hh"
#include "KalmanTests/inc/KalFitResult.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"

///using namespace CLHEP;

namespace mu2e {
  typedef std::vector<TrkStrawHit*>::iterator TSHI;

  FixedAmbigResolver::FixedAmbigResolver(fhicl::ParameterSet const& pset) : AmbigResolver(pset), 
     _neutralize(pset.get<bool>("Neutralize",true))
 {}

  FixedAmbigResolver::~FixedAmbigResolver() {}

  void
  FixedAmbigResolver::resolveTrk(KalFitResult& kfit) const {
// loop over all the hits
    TSHI ihit = kfit._hits.begin();
    while(ihit != kfit._hits.end()){
      TrkStrawHit* hit = *ihit++;
// set external error and don't allow the hit to auto-update its ambiguity
      hit->setExtErr     (AmbigResolver::_extErr);
      hit->setAmbigUpdate(false);
      if(_neutralize)hit->setAmbig(0);
    }
  }
}
