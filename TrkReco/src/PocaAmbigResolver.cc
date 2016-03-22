//
// class to resolve hit ambiguities one hit at a time, assuming a reasonable track
// fit as input
//
// $Id: PocaAmbigResolver.cc,v 1.2 2012/08/31 22:39:00 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/08/31 22:39:00 $
//
#include "TrkReco/inc/PocaAmbigResolver.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/KalmanTrack/KalSite.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include <vector>
#include <algorithm>
#include <functional>

namespace mu2e {
  typedef std::vector<TrkStrawHit*>::iterator TSHI;

  PocaAmbigResolver::PocaAmbigResolver(fhicl::ParameterSet const& pset, double ExtErr): 
    AmbigResolver(ExtErr)
  {}
  
  PocaAmbigResolver::~PocaAmbigResolver() {}
  
  bool 
  PocaAmbigResolver::resolveTrk(KalRep* krep) const {
    // init hit external errors
    initHitErrors(krep);
    // loop over all the hits
    TrkStrawHitVector tshv;
    convert(krep->hitVector(),tshv);
    for (auto ihit=tshv.begin();ihit!=tshv.end(); ++ihit){
      (*ihit)->setAmbigUpdate(true);
    }
    return true;
  }
}
