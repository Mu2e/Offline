//
// class to resolve hit ambiguities one hit at a time, assuming a reasonable track
// fit as input
//
// $Id: FixedAmbigResolver.cc,v 1.1 2012/05/14 19:20:02 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/05/14 19:20:02 $
//
#include "KalmanTests/inc/FixedAmbigResolver.hh"
#include "KalmanTests/inc/KalFit.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"

///using namespace CLHEP;

namespace mu2e {
  typedef std::vector<TrkStrawHit*>::iterator TSHI;

  FixedAmbigResolver::FixedAmbigResolver(fhicl::ParameterSet const& pset) : AmbigResolver(pset), 
     _neutralize(pset.get<bool>("Neutralize",false))
 {}

  FixedAmbigResolver::~FixedAmbigResolver() {}

  void
  FixedAmbigResolver::resolveTrk(TrkKalFit& kfit) const {
// loop over all the hits
    TSHI ihit = kfit._hits.begin();
    while(ihit != kfit._hits.end()){
      TrkStrawHit* hit = *ihit++;
      hit->setAmbigUpdate(false);
      if(_neutralize)hit->setAmbig(0);
    }
  }
}
