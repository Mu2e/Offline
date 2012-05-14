//
// class to resolve hit ambiguities one hit at a time, assuming a reasonable track
// fit as input
//
// $Id: PocaAmbigResolver.cc,v 1.1 2012/05/14 19:20:02 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/05/14 19:20:02 $
//
#include "KalmanTests/inc/PocaAmbigResolver.hh"
#include "KalmanTests/inc/KalFit.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTrack/KalRep.hh"
#include "KalmanTrack/KalSite.hh"
#include "KalmanTrack/KalHit.hh"
#include "TrkBase/TrkPoca.hh"
#include <vector>
#include <algorithm>
#include <functional>

///using namespace CLHEP;

namespace mu2e {
  typedef std::vector<TrkStrawHit*>::iterator TSHI;

  PocaAmbigResolver::PocaAmbigResolver(fhicl::ParameterSet const& pset) : AmbigResolver(pset)
  {}

  PocaAmbigResolver::~PocaAmbigResolver() {}

  void
  PocaAmbigResolver::resolveTrk(TrkKalFit& kfit) const {
// loop over all the hits
    TSHI ihit = kfit._hits.begin();
    while(ihit != kfit._hits.end()){
      TrkStrawHit* hit = *ihit++;
      hit->setAmbigUpdate(true);
    }
  }
}
