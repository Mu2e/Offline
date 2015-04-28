//
// Base class to resolve hit ambiguities
//
// $Id: AmbigResolver.cc,v 1.1 2012/05/14 19:20:02 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/05/14 19:20:02 $
//
#include "KalmanTests/inc/AmbigResolver.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTrack/KalRep.hh"
#include "KalmanTrack/KalSite.hh"
#include "KalmanTrack/KalHit.hh"
#include "KalmanTests/inc/KalFitResult.hh"
#include "TrkBase/TrkPoca.hh"
#include "difAlgebra/DifPoint.hh"
#include "difAlgebra/DifVector.hh"
#include <vector>
#include <algorithm>
#include <functional>

namespace mu2e {
  typedef std::vector<TrkStrawHit*>::iterator TSHI;
  typedef std::vector<TrkStrawHit*>::const_iterator TSHCI;
  typedef std::vector<KalSite*>::const_iterator KSI;

  AmbigResolver::AmbigResolver(fhicl::ParameterSet const& pset, double ExtErr, int Iter) {
    _extErr = ExtErr;
    _iter   = Iter;
  }

  AmbigResolver::~AmbigResolver() {}

//-----------------------------------------------------------------------------
// in the beginning of iteration set external hit errors to a constant
//-----------------------------------------------------------------------------
  void AmbigResolver::initHitErrors(KalFitResult& KRes) const {
    TrkStrawHit *hit;
    
    int nhits = KRes._hits.size();
    for (int i=0; i<nhits; ++i) {
      hit = KRes._hits.at(i);
      hit->setExtErr(_extErr);
    }
  }

  const TrkSimpTraj*
  AmbigResolver::findTraj(std::vector<TrkStrawHit*> const& phits, const KalRep* krep) const {
    const TrkSimpTraj* retval(0);
// if the fit is valid, use the full fit result.  Otherwise, use the reference traj
    if(!krep->fitValid()){
      double locdist;
      retval = krep->referenceTraj()->localTrajectory(phits[0]->fltLen(),locdist);
    } else {
// find the range of these hits in the KalRep site vector
      std::vector<KalSite*> const& sites = krep->siteList();
      KSI first = sites.end();
      KSI last = sites.begin();
      for(TSHCI ihit = phits.begin();ihit != phits.end();++ihit){
	if((*ihit)->isActive()){
	  const KalHit* hitsite = krep->findHotSite(*ihit);
	  // find the index to this site
	  KSI ifnd = std::find(sites.begin(),sites.end(),hitsite);
	  if(ifnd != sites.end()) {
	    if(ifnd < first){
	      first = ifnd;
	    }
	    if(ifnd > last){
	      last = ifnd;
	    }
	  }
	}
      }
// create a trajectory from the fit which excludes this set of hits
      static TrkSimpTraj* straj = krep->seed()->clone();
      if(krep->smoothedTraj(first,last,straj)){
	retval = straj;
      } else {
	double locdist;
	retval = krep->referenceTraj()->localTrajectory(phits[0]->fltLen(),locdist);
      }
    }
    return retval;
  }


}
