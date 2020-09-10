// class to resolve hit ambiguities by panel, assuming a reasonable track
// fit as input
//
//
#include "TrkReco/inc/PanelAmbigResolver.hh"
#include "TrkReco/inc/PanelStateIterator.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/TrkBase/TrkT0.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/KalmanTrack/KalSite.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/difAlgebra/DifPoint.hh"
#include "BTrk/difAlgebra/DifVector.hh"
#include <vector>
#include <algorithm>
#include <functional>
#include <iostream>
// art
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
using CLHEP::Hep3Vector;
using CLHEP::HepSymMatrix;
using CLHEP::HepVector;
namespace mu2e {
  namespace PanelAmbig {
    // functor for sorting by panel.  Note that PanelId uniquely defines a panel
    struct panelcomp : public std::binary_function<TrkStrawHit*, TrkStrawHit*, bool> {
      bool operator()(TrkStrawHit* x, TrkStrawHit* y) { return x->straw().id().getPanelId() < y->straw().id().getPanelId(); }
    };

    // functor to sort panel results by chisquared
    struct resultcomp : public std::binary_function<PanelResult const&, PanelResult const&,bool> {
      bool operator()(PanelResult const& a, PanelResult const& b) { return a._chisq < b._chisq; }
    };

    typedef TrkStrawHitVector::iterator TSHI;
    typedef TrkStrawHitVector::const_iterator TSHCI;


    PanelAmbigResolver::PanelAmbigResolver(fhicl::ParameterSet const& pset, double tmpErr, size_t iter): 
      AmbigResolver(tmpErr),
      _minsep(pset.get<double>("minChisqSep",4.0)),
      _inactivepenalty(pset.get<double>("inactivePenalty",16.0)),
      _nullpenalty(pset.get<double>("NullHitPenalty",0.0)),
      _penaltyres(pset.get<double>("PenaltyResolution",0.25)),
      _trkpenaltyres(pset.get<double>("TrackPenaltyResolution",0.5)),
      _addtrkpos(pset.get<bool>("AddTrackPositionConstraint",true)),
      _maxhitu(pset.get<double>("MaximumHitU",8.0)),
      _fixunallowed(pset.get<bool>("FixUnallowedHitStates",true)),
      _maxnpanel(pset.get<unsigned>("MaxHitsPerPanel",8)),
      _diag(pset.get<int>("DiagLevel",0))
    {
      double nullerr = pset.get<double>("ExtraNullAmbigError",0.0);
      _nullerr2 = nullerr*nullerr;
      std::vector<int> allowed = pset.get< std::vector<int> >("AllowedHitStates");
      for(std::vector<int>::iterator ial = allowed.begin();ial != allowed.end();++ial){
	if(*ial >= HitState::negambig && *ial <= HitState::inactive){
	  _allowed.push_back(HitState(static_cast<HitState::TSHState>(*ial)));
	} else
	  throw cet::exception("RECO")<<"mu2e::PanelAmbigResolver: illegal state" << std::endl;
      }
      // if requested, setup diagnostics
      if(_diag > 0){
	art::ServiceHandle<art::TFileService> tfs;
	char title[40];
	snprintf(title,40,"padiag_%lu",iter);
	_padiag=tfs->make<TTree>(title,"Panel Ambiguity Resolution Diagnostics");
	_padiag->Branch("nhits",&_nrhits,"nhits/I");
	_padiag->Branch("nactive",&_nactive,"nactive/I");
	_padiag->Branch("nres",&_nres,"nres/I");
	_padiag->Branch("results",&_results);

	snprintf(title,40,"pudiag_%lu",iter);
	_pudiag=tfs->make<TTree>(title,"Panel U position Diagnostics");
	_pudiag->Branch("nhits",&_nuhits,"nhits/I");
	_pudiag->Branch("tupos",&_tupos,"tupos/F");
	_pudiag->Branch("tuerr",&_tuerr,"tuerr/F");
	_pudiag->Branch("uinfo",&_uinfo);
      }
    }

    PanelAmbigResolver::~PanelAmbigResolver() {}

    bool PanelAmbigResolver::resolveTrk(KalRep* krep) const {
      bool retval(false); // assume nothing changes
      // initialize penalty errors
      initHitErrors(krep);
      // sort by panel
      TrkStrawHitVector tshv;
      convert(krep->hitVector(),tshv);
      std::sort(tshv.begin(),tshv.end(),panelcomp());
      // collect hits in the same panel
      auto ihit=tshv.begin();
      while(ihit!=tshv.end()){
	PanelId pid = (*ihit)->straw().id().getPanelId();
	(*ihit)->setTemperature(AmbigResolver::_tmpErr);
	TrkStrawHitVector phits;
	auto jhit=ihit;
	while(jhit != tshv.end() && (*jhit)->straw().id().getPanelId() == pid){
	  phits.push_back(*jhit++);
	}
	// resolve the panel hits
	retval |= resolvePanel(phits,krep);
	ihit = jhit;
      }
      return retval;
    }

    bool PanelAmbigResolver::resolvePanel(TrkStrawHitVector& phits,KalRep* krep) const {
      bool retval(false); // assume nothing changes
      // sort hits for this panel
      std::sort(phits.begin(),phits.end(),hitsort());
      // fill panel information
      PanelInfo pinfo;
      if(fillPanelInfo(phits,krep,pinfo)){
	// loop over all ambiguity/activity states for this panel
	PanelStateIterator psi(pinfo._uinfo,_allowed);
	PRV results;
	do {
	  // for each state, fill the result of the 1-dimensional optimization
	  PanelResult result(psi.current());
	  fillResult(pinfo,krep->t0(),result);
	  if(result._status == 0)results.push_back(result);
	} while(psi.increment());
	if(results.size() > 0){
	  // sort the results to have lowest chisquard first
	  std::sort(results.begin(),results.end(),resultcomp());
	  // for now, set the hit state according to the best result.  In future, maybe we want to treat
	  // cases with different ambiguities differently from inactive hits
	  retval |= setHitStates(results[0]._state,phits);
	  // if the chisq difference between patterns is negligible, inflate the errors of the
	  // hit which changes
	  size_t nhits = results[0]._state.size();
	  size_t ires(1);
	  while (ires < results.size() && results[ires]._chisq - results[0]._chisq < _minsep){
	    for(size_t ihit=0;ihit<nhits;++ihit){
	      if(results[ires]._state[ihit] != results[0]._state[ihit]){
		phits[ihit]->setPenalty(_penaltyres);
	      }
	    }
	    ++ires;
	  }
	}
	if( _diag > 1 ) {
	  _nuhits = _nrhits = pinfo._uinfo.size();
	  _nactive = 0;
	  for(auto const& ishi : pinfo._uinfo) {
	    if(ishi._active)++_nactive;
	  }
	  _nres = results.size();
	  _results = results;
	  _padiag->Fill();
	  //
	  _tupos = pinfo._tupos;
	  _tuerr = pinfo._tuerr;
	  _uinfo = pinfo._uinfo;
	  _pudiag->Fill();
	}
      } else 
	std::cout << "PanelAmbigResolver: Panel with " << phits.size() << " hits has no usable info" << std::endl;

      return retval;
    }

    bool PanelAmbigResolver::fillPanelInfo(TrkStrawHitVector const& phits, const KalRep* krep, PanelInfo& pinfo) const {
      bool retval(false);
      // find the best trajectory we can local to these hits, but excluding their information ( if possible).
      const TrkSimpTraj* straj = findTraj(phits,krep);
      if(straj != 0){
	// find a reference point on this traj using POCA to the first good hit
	double gdist;
	Hep3Vector wdir; 
	HepPoint wpos;
	for(const TrkStrawHit* refhit : phits) {
	  TrkPoca tpoca(*straj,refhit->fltLen(),*refhit->hitTraj(),refhit->hitLen());
	  if(tpoca.status().success() ) { //&& abs(tpoca.doca()) < _maxhitu){
	    retval = true;
	    gdist = tpoca.flt1();
	    wdir = refhit->hitTraj()->direction(tpoca.flt2()); 
	    wpos = refhit->hitTraj()->position(tpoca.flt2()); 
	    break;
	  }
	}
	if(retval){
	  // find the trajectory position and direction information at this distance
	  DifPoint tpos;
	  DifVector tdir;
	  straj->getDFInfo2(gdist, tpos, tdir);
	  // straw information; cast to dif even though the derivatives are 0
	  Hep3Vector wposv(wpos.x(),wpos.y(),wpos.z());
	  DifVector delta = DifVector(wposv) - tpos;
	  // compute the U direction (along the measurement, perp to the track and wire
	  // The sign is chosen to give positive and negative ambiguity per BaBar convention
	  DifVector dudir = cross(DifVector(wdir),tdir).normalize();
	  Hep3Vector udir(dudir.x.number(),dudir.y.number(),dudir.z.number());
	  // compute the track constraint on u, and it's error.  The dif algebra part propagates the error
	  DifNumber trku = dudir.dot(delta);
	  pinfo._udir = udir;
	  Hep3Vector tposv(tpos.x.number(),tpos.y.number(),tpos.z.number());
	  pinfo._tupos = udir.dot(tposv);
	  pinfo._tuerr = trku.error();
	  // degrade the track information
	  pinfo._tuwt = 1.0/(pinfo._tuerr*pinfo._tuerr + _trkpenaltyres*_trkpenaltyres);
	  // now, project the hits onto u, WRT the projected track position. 
	  for(auto ihit = phits.begin();ihit != phits.end();++ihit){
	    TSHUInfo uinfo(*ihit,udir,tpos.hepPoint());
	    // mask off hits with wire u values outside the limits, and (optionally) those whose
	    // initial state isn't allowed
	    if(fabs(uinfo._upos) > _maxhitu)
	      uinfo._use = TSHUInfo::unused;
	    if(_fixunallowed){
	      auto ifnd = std::find(_allowed.begin(),_allowed.end(),uinfo._hstate);
	      if(ifnd == _allowed.end())
		uinfo._use = TSHUInfo::fixed;
	    }
	    // limit the # of hits considered
	    if(pinfo._nfree >= _maxnpanel){
	      std::cout << "PanelAmbigResolver: Maximum number of hits/panel exceeded: truncating" << std::endl;
	      uinfo._use = TSHUInfo::unused;
	    }
	    if(uinfo._use == TSHUInfo::free)++pinfo._nfree;
	    if(uinfo._use != TSHUInfo::unused)++pinfo._nused;
	    pinfo._uinfo.push_back(uinfo);
	  }
	}
      } else {
	throw cet::exception("RECO")<<"mu2e::PanelAmbigResolver: no trajectory" << std::endl;
      }
      return retval;
    }

    void PanelAmbigResolver::fillResult(PanelInfo const& pinfo,TrkT0 const& t0, PanelResult& result) const {
      // initialize the sums
      double wsum(0.0);
      double uwsum(0.0);
      double vwsum(0.0);
      double uuwsum(0.0);
      double vvwsum(0.0);
      double uvwsum(0.0);
      double chi2penalty(0.0);
      // loop over the straw hit info and accumulate the sums used to compute chisquared
      size_t ntsh = pinfo._uinfo.size();
      // consistency check
      if(ntsh != result._state.size()) 
	throw cet::exception("RECO")<<"mu2e::PanelAmbigResolver: inconsistent hits" << std::endl;
      for(size_t itsh=0;itsh<ntsh;++itsh){
	HitState const& tshs = result._state[itsh];
	TSHUInfo const& tshui = pinfo._uinfo[itsh];
	// compute chisquared for selected hits
	if(tshui._use == TSHUInfo::free || tshui._use == TSHUInfo::fixed){
	  // compare this state to the original, record any differences
	  if(tshs != tshui._hstate)
	    result._statechange |= (itsh << 1);
	  // compute u 
	  if(tshs._state != HitState::inactive){
	    double w = tshui._uwt;
	    double r = tshui._dr;
	    double v = tshui._dv;
	    // sign for ambiguity
	    if(tshs._state == HitState::negambig){
	      r *= -1;
	      v *= -1;
	    } else if(tshs._state == HitState::noambig){
	      r = 0.; // inactive hits don't depend on time
	      v = 0.;
	      w = 1.0/(1.0/w + _nullerr2); // increase the error on 0 ambiguity hits
	      chi2penalty += _nullpenalty;
	    }
	    double u = tshui._upos + r;
	    wsum += w;
	    uwsum += u*w;
	    vwsum += v*w;
	    uuwsum += u*u*w;
	    vvwsum += v*v*w;
	    uvwsum += u*v*w;
	  } else // penalize inactive hits
	    chi2penalty += _inactivepenalty;
	}
      }
      if(pinfo._nused > 0){
	// propogate t0 uncertainty.  This is a constraint centered at the
	// current value of t0, unit derivative
	double t0wt = 1.0/(t0._t0err*t0._t0err);
	vvwsum += t0wt; // t0 has unit derrivative
	// optionally add track position constraint if there are multiple hits.  It is like a hit, but with r=v=0
	// NB, since the track position is defined to be 0, it adds only to the weight
	if(_addtrkpos || pinfo._nused == 1)wsum += pinfo._tuwt;
	// now compute the scalar, vector, and tensor terms for the 2nd-order chisquared expansion
	double alpha = uuwsum;
	HepVector beta(2);
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
#ifdef DEBUG
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
#endif
	  // add a penalty term (if any) for this particular pattern.  Still to be written, FIXME!!!
	  //  addPatternPenalty();
	}
      } else {
      // if there are no active hits, flag the status
	result._status = -100;
      }
    }

  } // PanelAmbig namespace
} // mu2e namespace
