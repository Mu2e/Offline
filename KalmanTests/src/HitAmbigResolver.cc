//
// class to resolve hit ambiguities one hit at a time, assuming a reasonable track
// fit as input
//
// $Id: HitAmbigResolver.cc,v 1.3 2012/08/31 22:39:00 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/08/31 22:39:00 $
//
#include "KalmanTests/inc/HitAmbigResolver.hh"
#include "KalmanTests/inc/KalFitResult.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
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

  HitAmbigResolver::HitAmbigResolver(fhicl::ParameterSet const& pset, double ExtErr, int Iter) : 
    AmbigResolver(pset,ExtErr,Iter),
    _mindrift(pset.get<double>("HitMinDrift",0.1)),
    _zeropenalty(pset.get<double>("ZeroDriftPenalty",0.05)),
    _penalty(pset.get<bool>("HitAmbigPenalty",false)),
    _expnorm(pset.get<double>("HitExpNorm",0.03907)),
    _lambda(pset.get<double>("HitLambda",0.1254)),
    _offset(pset.get<double>("HitOffset",0.073)),
    _slope(pset.get<double>("HitSlope",-0.002374))
  {}

  HitAmbigResolver::~HitAmbigResolver() {}

  void
  HitAmbigResolver::resolveTrk(KalFitResult& kfit, int Final) const {

    initHitErrors(kfit);

// loop over all the hits
    TSHI ihit = kfit._hits.begin();
    while(ihit != kfit._hits.end()){
      TrkStrawHit* hit = *ihit++;
// don't allow the hit to auto-update its ambiguity
      hit->setAmbigUpdate(false);
// get the drift radius
      double rdrift = hit->driftRadius();
      if(rdrift <= _mindrift){
	hit->setAmbig(0);
	hit->setPenalty(_zeropenalty);
      } else {
// find the best trajectory we can local to these hits, but excluding their information ( if possible).
	std::vector<TrkStrawHit*> hits;
	hits.push_back(hit);
	const TrkDifTraj* traj = findTraj(hits,kfit._krep);
	// compute POCA to this traj
	TrkPoca poca(*traj,hit->fltLen(),*hit->hitTraj(),hit->hitLen());
	if(poca.status().success()){
	  // set the ambiguity if allowed, based on the sign of DOCA
	  int newamb = poca.doca() > 0 ? 1 : -1;
	  hit->setAmbig(newamb);
	  // based on the drift distance, set the penalty error based on the a-priori function.
	  if(_penalty){
	    double perr = penaltyError(hit->driftRadius());
	    hit->setPenalty(perr);
	  }
	}
      }
    }
  }

  double
  HitAmbigResolver::penaltyError(double rdrift) const {
    // model is of an exponential plus a linear.
    double frac = _expnorm*exp(-rdrift/_lambda)/_lambda + _offset + _slope*rdrift;
    // the penalty term for a discrete ambiguity error depends only on the mis-assignment probability and the drift distance
    static double sqrt2 = sqrt(2.0);
    double perr = sqrt2*rdrift*sqrt(frac);
    return perr;
  }
}
