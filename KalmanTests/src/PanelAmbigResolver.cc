//
// class to resolve hit ambiguities by panel, assuming a reasonable track
// fit as input
//
// $Id: PanelAmbigResolver.cc,v 1.1 2012/05/14 19:20:02 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/05/14 19:20:02 $
//
#include "KalmanTests/inc/PanelAmbigResolver.hh"
#include "KalmanTests/inc/KalFit.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTrack/KalRep.hh"
#include "KalmanTrack/KalSite.hh"
#include "KalmanTrack/KalHit.hh"
#include "TrkBase/TrkPoca.hh"
#include "difAlgebra/DifPoint.hh"
#include "difAlgebra/DifVector.hh"
#include <vector>
#include <algorithm>
#include <functional>

///using namespace CLHEP;

namespace mu2e {
// functor for sorting by panel.  Note that SectorId uniquely defines a panel
  struct panelcomp : public std::binary_function<TrkStrawHit*, TrkStrawHit*, bool> {
    bool operator()(TrkStrawHit* x, TrkStrawHit* y) { return x->straw().id().getSectorId() < y->straw().id().getSectorId(); }
  };

// functor to sort panel results by chisquared
  struct resultcomp : public std::binary_function<PanelResult const&, PanelResult const&,bool> {
    bool operator()(PanelResult const& a, PanelResult const& b) { return a._chisq < b._chisq; }
  };

  TSHUInfo::TSHUInfo(const TrkStrawHit* tsh,CLHEP::Hep3Vector const& udir, CLHEP::Hep3Vector const& uorigin): _tsh(tsh) {
    CLHEP::Hep3Vector dstraw = tsh->straw().getMidPoint()-uorigin;
    _upos = udir.dot(dstraw);
    // note that the t0 component of the error scales coherently between the hits, so we use only the intrinsic error
    double uerr = tsh->hitErr();
    _uwt = 1.0/(uerr*uerr);
  }

  typedef std::vector<TrkStrawHit*>::iterator TSHI;
  typedef std::vector<TrkStrawHit*>::const_iterator TSHCI;
  typedef std::vector<KalSite*>::const_iterator KSI;

  PanelAmbigResolver::PanelAmbigResolver(fhicl::ParameterSet const& pset) : AmbigResolver(pset),
    _minsep(pset.get<double>("minChisqSep",6.0)),
    _inactivepenalty(pset.get<double>("inactivePenalty",16.0)),
    _penaltyres(pset.get<double>("PenaltyResolution",0.5)),
    _setactivity(pset.get<bool>("sethitactivity",true))
  {
    std::vector<int> allowed = pset.get< std::vector<int> >("AllowedHitStates");
    for(std::vector<int>::iterator ial = allowed.begin();ial != allowed.end();++ial){
      if(*ial >= TrkStrawHitState::noambig && *ial < TrkStrawHitState::nstates){
	_allowed.push_back(static_cast<TrkStrawHitState::TSHState>(*ial));
      } else
	throw cet::exception("RECO")<<"mu2e::PanelAmbigResolver: illegal state" << std::endl;
    }
  }

  PanelAmbigResolver::~PanelAmbigResolver() {}

  void 
  PanelAmbigResolver::resolveTrk(TrkKalFit& kfit) const {
// sort by panel
    std::sort(kfit._hits.begin(),kfit._hits.end(),panelcomp());
// collect hits in the same panel
    TSHI ihit = kfit._hits.begin();
    while(ihit != kfit._hits.end()){ 
      SectorId pid = (*ihit)->straw().id().getSectorId();
      TSHI jhit = ihit;
      std::vector<TrkStrawHit*> phits;
      while(jhit != kfit._hits.end() && (*jhit)->straw().id().getSectorId() == pid){
	phits.push_back(*jhit++);
      }
// resolve the panel hits
      resolvePanel(phits,kfit);
      ihit = jhit;
    }
  }

  void
  PanelAmbigResolver::resolvePanel(std::vector<TrkStrawHit*>& phits,TrkKalFit& kfit) const {
// fill panel information
    PanelInfo pinfo;
    fillPanelInfo(phits,kfit._krep,pinfo);
// loop over ambiguity/activity states for these hits
    TrkStrawHitStateVector tshsv(phits,_allowed);
    tshsv.reset();
    std::vector<PanelResult> results;
    do {
// for each state, fill the result of the 1-dimensional optimization
      PanelResult result(tshsv);
      fillResult(pinfo,kfit._t0,result);
      if(result._status == 0)results.push_back(result);
    } while(tshsv.increment());
// sort the results to have lowest chisquard first
    std::sort(results.begin(),results.end(),resultcomp());
// compare the chisquared of the top 2 values: if these are separated by more than the minimum,
// we consider this panel fully resolved
    if(results.size() == 1 || (results.size() > 1 && results[1]._chisq - results[0]._chisq > _minsep) ){
      results[0]._state.setHitStates(_setactivity);
// apply the state of the solution
    } else if(results.size()>1){
// for now, set the hit state according to the best result.  In future, maybe we want to treat
// cases with different ambiguities differently from inactive hits
      results[0]._state.setHitStates(_setactivity);
// determine the state of hits one-by-one: if they are different between patterns consistent within the
// minimum, inflate their errors
      size_t nhits = results[0]._state.nHits();
      unsigned ires(1);
      do {
	for(size_t ihit=0;ihit<nhits;++ihit){
	  if(results[ires]._state.states()[ihit] != results[0]._state.states()[ihit]){
	    results[ires]._state.states()[ihit].hit()->setPenalty(_penaltyres);
	  }
	}
	++ires;
      } while (ires < results.size() && results[ires]._chisq - results[0]._chisq < _minsep);
    }
  }

  bool
  PanelAmbigResolver::fillPanelInfo(std::vector<TrkStrawHit*> const& phits, const KalRep* krep, PanelInfo& pinfo) const {
    bool retval(false);
// find the best trajectory we can local to these hits, but excluding their information ( if possible).
    const TrkSimpTraj* straj = findTraj(phits,krep);
    if(straj != 0){
// find POCA to the first hit
      const TrkStrawHit* firsthit = phits[0];
      TrkPoca tpoca(*straj,firsthit->fltLen(),*firsthit->hitTraj(),firsthit->hitLen());
      if(tpoca.status().success()){
	retval = true;
	// find the position and direction information at the middle of this range.
	DifPoint tpos;
	DifVector tdir;
	straj->getDFInfo2(tpoca.flt1(), tpos, tdir);
	// straw information; cast to dif even though the derivatives are 0
	CLHEP::Hep3Vector wdir = firsthit->straw().getDirection();
	CLHEP::Hep3Vector wpos = firsthit->straw().getMidPoint();
	DifVector delta = DifVector(wpos) - tpos;
	// compute the U direction (along the measurement, perp to the track and wire
	DifVector dudir = cross(tdir,DifVector(wdir));
	CLHEP::Hep3Vector udir(dudir.x.number(),dudir.y.number(),dudir.z.number());
	// compute the track constraint on u, and it's error.  The dif algebra part propagates the error
	DifNumber trku = dudir.dot(delta);
	pinfo._utpos = trku.number();
	double uterr = trku.error();
	pinfo._utwt = 1.0/(uterr*uterr);
	// now, compute the projections of all the hits (including inactive) WRT the firsthit hit.  This creates an effective 1-d measure for all hits
	for(TSHCI ihit = phits.begin();ihit != phits.end();++ihit){
	  pinfo._uinfo.push_back(TSHUInfo(*ihit,udir,wpos));
	}
      }
    } else {
      throw cet::exception("RECO")<<"mu2e::PanelAmbigResolver: no trajectory" << std::endl;
    }
    return retval;
  }

  void
  PanelAmbigResolver::fillResult(PanelInfo const& pinfo,TrkT0 const& t0, PanelResult& result) const {
// initialize the sums
    double wsum(0.0);
    double uwsum(0.0);
    double vwsum(0.0);
    double uuwsum(0.0);
    double vvwsum(0.0);
    double uvwsum(0.0);
    double chi2penalty(0.0);
// loop over the straw hit info and state in parallel, and accumulate the sums
    size_t ntsh = pinfo._uinfo.size();
    std::vector<TrkStrawHitState> const& hitstates = result._state.states();
    if(ntsh != hitstates.size())  
      throw cet::exception("RECO")<<"mu2e::PanelAmbigResolver: inconsistent hits" << std::endl;
    for(size_t itsh=0;itsh<ntsh;++itsh){
// only accumulate if the state is active
      TrkStrawHitState const& tshs = hitstates[itsh];
      TSHUInfo const& tshui = pinfo._uinfo[itsh];
// consistency check
      if(tshs.hit() != tshui._tsh)
	throw cet::exception("RECO")<<"mu2e::PanelAmbigResolver: inconsistent hits" << std::endl;
      if(tshs.state() != TrkStrawHitState::inactive){
	double w = tshui._uwt;
	double r = tshui._tsh->driftRadius();
	double v = tshui._tsh->driftVelocity();
// sign for ambiguity
	if(tshs.state() == TrkStrawHitState::negambig){
	  r *= -1;
	  v *= -1;
	}
	double u = tshui._upos + r;
	wsum += w;
	uwsum += u*w;
	vwsum += v*w;
	uuwsum += u*u*w;
	vvwsum += v*v*w;
	uvwsum += u*v*w;
      } else
	chi2penalty += _inactivepenalty;
    }
// add track consistency information; position and t0
    double t0wt = 1.0/(t0._t0err*t0._t0err);
    vvwsum += t0wt; // t0 has unit derrivative
// track position is like a hit, but with r=v=0
    wsum += pinfo._utwt;
    uwsum += pinfo._utpos*pinfo._utwt;
    uuwsum += pinfo._utpos*pinfo._utpos*pinfo._utwt;
// now compute the scalar, vector, and tensor term for the 2nd-order chisquared expansion
    double alpha = uuwsum;
    CLHEP::HepVector beta(2);
    beta(1) = uwsum;
    beta(2) = uvwsum;
    HepSymMatrix gamma(2);
    gamma.fast(1,1) = wsum;
    gamma.fast(2,2) = vvwsum;
    gamma.fast(1,2) = vwsum;
// invert and solve
    gamma.invert(result._status);
    if(result._status == 0){
      result._dcov = gamma;
      result._delta = gamma * beta;
      result._chisq = alpha - gamma.similarity(beta) + chi2penalty;
// debug printout
      double g11 = gamma.fast(1,1);
      double g22 = gamma.fast(2,2);
      double g12 = gamma.fast(1,2);
      double g21 = gamma.fast(2,1);
      double b1 = beta(1);
      double b2 = beta(2);
      double d1 = result._delta(1);
      double d2 = result._delta(2);
      double sim = gamma.similarity(beta);
      double test = alpha - sim;
      double test2 = alpha -2*dot(beta,result._delta) + gamma.similarity(result._delta);
// add a penalty term (if any) for this particular pattern.  Still to be written, FIXME!!!
//  addPatternPenalty();
    }
  }

}
