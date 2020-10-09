//
// Base class to resolve hit ambiguities
//
//
#include "TrkReco/inc/AmbigResolver.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/KalmanTrack/KalSite.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/difAlgebra/DifPoint.hh"
#include "BTrk/difAlgebra/DifVector.hh"
#include <vector>
#include <algorithm>
#include <functional>

namespace mu2e {
  typedef std::vector<KalSite*>::const_iterator KSI;

  AmbigResolver::AmbigResolver(double tmpErr) : _tmpErr(tmpErr) {}

  AmbigResolver::~AmbigResolver() {}

//-----------------------------------------------------------------------------
// in the beginning of iteration set external hit errors to a constant
//-----------------------------------------------------------------------------
  void AmbigResolver::initHitErrors(KalRep* krep) const {
// get hits and cast to TrkStrawHits
    TrkStrawHitVector tshv;
    convert(krep->hitVector(),tshv);
    for (auto itsh=tshv.begin();itsh!=tshv.end(); ++itsh){
      (*itsh)->setTemperature(_tmpErr);
    }
  }

  const TrkSimpTraj*
  AmbigResolver::findTraj(std::vector<TrkStrawHit*> const& phits, const KalRep* krep) const {
    const TrkSimpTraj* retval(0);
// if the fit is valid, use the full fit result.
    if(krep->fitValid()){
// find the range of these hits in the KalRep site vector
      std::vector<KalSite*> const& sites = krep->siteList();
      KSI first = sites.begin();
      while(first != sites.end() && (*first)->globalLength() < phits.front()->fltLen() ) {
	++first;
      }
      KSI last = first;
      while(last != sites.end() && (*last)->globalLength() < phits.back()->fltLen() ){
	++last;
      }
      // back off one
      if(first != sites.begin())--first;
      if(last == sites.end())--last;
// create a trajectory from the fit which excludes this set of hits
// Use of static is memory-efficient but not threadsafe FIXME!!!
      static TrkSimpTraj* straj = krep->seed()->clone();
      if(krep->smoothedTraj(first,last,straj)){
	retval = straj;
      } 
    }
//  Otherwise, use the reference traj at the center of these hits
    if(retval == 0){
      double locdist;
      double gdist = 0.5*(phits.front()->fltLen()+phits.back()->fltLen());
      retval = krep->referenceTraj()->localTrajectory(gdist,locdist);
    }
    return retval;
  }


}
