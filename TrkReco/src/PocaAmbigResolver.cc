//
// class to resolve hit ambiguities one hit at a time, assuming a reasonable track
// fit as input
//
// $Id: PocaAmbigResolver.cc,v 1.2 2012/08/31 22:39:00 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/08/31 22:39:00 $
//
#include "TrkReco/inc/PocaAmbigResolver.hh"
#include "TrkReco/inc/KalFitResult.hh"
#include "Mu2eBTrk/inc/TrkStrawHit.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/KalmanTrack/KalSite.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include <vector>
#include <algorithm>
#include <functional>

///using namespace CLHEP;

namespace mu2e {
  typedef std::vector<TrkStrawHit*>::iterator TSHI;

  PocaAmbigResolver::PocaAmbigResolver(fhicl::ParameterSet const& pset, double ExtErr): 
    AmbigResolver(ExtErr)
  {}
  
  PocaAmbigResolver::~PocaAmbigResolver() {}
  
  void
  PocaAmbigResolver::resolveTrk(KalFitResult& kfit) const {
    // init hit external errors
    initHitErrors(kfit);
    // loop over all the hits
    TSHI ihit = kfit._hits.begin();
    while(ihit != kfit._hits.end()){
      TrkStrawHit* hit = *ihit++;
      hit->setAmbigUpdate(true);
    }
  }
}
