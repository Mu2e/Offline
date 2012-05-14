//
// Class to perform BaBar Kalman fit
//
// $Id: KalFit.cc,v 1.26 2012/05/14 19:20:02 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/05/14 19:20:02 $
//

// the following has to come before other BaBar includes
#include "BaBar/BaBar.hh"
#include "KalmanTests/inc/KalFit.hh"
#include "KalmanTests/inc/PanelAmbigResolver.hh"
#include "KalmanTests/inc/PocaAmbigResolver.hh"
#include "KalmanTests/inc/HitAmbigResolver.hh"
#include "KalmanTests/inc/FixedAmbigResolver.hh"
//geometry
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "BFieldGeom/inc/BFieldConfig.hh"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
// tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
// BaBar
#include "KalmanTrack/KalHit.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/HelixTraj.hh"
#include "TrkBase/TrkHotListFull.hh"
#include "TrkBase/TrkHelixUtils.hh"
#include "TrkBase/TrkMomCalculator.hh"
#include "TrkBase/TrkPoca.hh"
#include "BaBar/ErrLog.hh"
#include "BField/BFieldFixed.hh"
#include "DetectorModel/DetIntersection.hh"
#include "DetectorModel/DetMaterial.hh"
#include "BField/BFieldFixed.hh"
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <memory>


using namespace std; 

namespace mu2e 
{
// convert speed of light in mm/nsec
  const double KalFit::_vlight = CLHEP::c_light;
// comparison functor for ordering hits
  struct fltlencomp : public binary_function<TrkStrawHit*, TrkStrawHit*, bool> {
          bool operator()(TrkStrawHit* x, TrkStrawHit* y) { return x->fltLen() < y->fltLen(); }
  };

  void TrkKalFit::deleteTrack() {

    if(_trk != 0 && _krep != 0){
      _hits.clear(); delete _trk; _trk=0; _krep = 0; 
//      std::cout << "deleting fit with track " << std::endl;
    } else {
// if there's no fit, we need to delete the hits
//      std::cout << "deleting " << _hits.size() << " hits from fit without track " << std::endl;
      std::for_each(_hits.begin(),_hits.end(),babar::Collection::DeleteObject());
    }
  } 


// construct from a parameter set  
  KalFit::KalFit(fhicl::ParameterSet const& pset) : _bfield(0),
    _debug(pset.get<int>("debugLevel",0)),
    _fieldcorr(pset.get<bool>("fieldCorrections",false)),
    _material(pset.get<bool>("material",true)),
//    _ambigflip(pset.get<bool>("ambigflip",false)),
    _weedhits(pset.get<bool>("weedhits",true)),
    _updatet0(pset.get<bool>("updateT0",true)),
    _removefailed(pset.get<bool>("RemoveFailedFits",true)),
    _t0tol(pset.get<double>("t0Tolerance",1.0)),
    _maxhitchi(pset.get<double>("maxhitchi",5.0)),
    _maxiter(pset.get<unsigned>("maxiter",10)),
    _mingap(pset.get<double>("mingap",0.1)),
    _minnstraws(pset.get<unsigned>("minnstraws",15)),
    _minndof(pset.get<unsigned>("minNDOF",10)),
    _maxweed(pset.get<unsigned>("maxweed",10)),
    _herr(pset.get< vector<double> >("hiterr")),
    _ssmear(pset.get<double>("seedsmear",1e6)),
    _t0errfac(pset.get<double>("t0ErrorFactor",1.2)),
    _mint0doca(pset.get<double>("minT0DOCA",-0.2)),
    _maxdriftpull(pset.get<double>("maxDriftPull",10)),
    _t0nsig(pset.get<double>("t0window",2.5)),
    _fitpart(pset.get<int>("fitparticle",PdtPid::electron)),
    _t0strategy((t0Strategy)pset.get<int>("t0Strategy",median)),
    _ambigstrategy(pset.get< vector<int> >("ambiguityStrategy"))
  {
//      _herr.push_back( pset.get< double >("hiterr") );
  // make sure we have at least one entry for additional errors
      if(_herr.size() <= 0) throw cet::exception("RECO")<<"mu2e::KalFit: no hit errors specified" << endl;
      if(_herr.size() != _ambigstrategy.size()) throw cet::exception("RECO")<<"mu2e::KalFit: inconsistent ambiguity resolution" << endl;
      _kalcon = new KalContext;
      _kalcon->setBendSites(_fieldcorr);
      _kalcon->setMaterialSites(_material);
    // depending on the ambiguity strategy, either allow or not the hit to change it itself
      _kalcon->setMaxIterations(_maxiter);
      _kalcon->setMinGap(_mingap); // minimum separation between sites when creating trajectory pieces
      // these are currently fixed, they should be set as parameters and re-optimized FIXME!!!!
      _kalcon->setMaxIntersections(0);
      _kalcon->setMaxDMom(10);
      _kalcon->setSmearFactor(_ssmear);
      _kalcon->setMinDOF(_minndof,TrkEnums::bothView);
      _kalcon->setMinDOF(_minndof,TrkEnums::xyView);
      _kalcon->setMinDOF(0,TrkEnums::zView);
      _kalcon->setIntersectionTolerance(100);
      _kalcon->setMaxMomDiff(1.0); // 1 MeV
      _kalcon->setTrajBuffer(0.01); // 10um
      _kalcon->setDefaultType((PdtPid::PidType)_fitpart);
// construct the ambiguity resolvers
      for(size_t iambig=0;iambig<_ambigstrategy.size();++iambig){
	switch (_ambigstrategy[iambig] ){
	  case fixedambig: default:
	    _ambigresolver.push_back(new FixedAmbigResolver(pset));
	    break;
	  case hitambig:
	    _ambigresolver.push_back(new HitAmbigResolver(pset));
	    break;
	  case panelambig:
	    _ambigresolver.push_back(new PanelAmbigResolver(pset));
	    break;
	  case pocaambig:
	    _ambigresolver.push_back(new PocaAmbigResolver(pset));
	    break;
	}
      }

    }

  KalFit::~KalFit(){
    delete _kalcon;
    for(size_t iambig=0;iambig<_ambigresolver.size();++iambig){
      delete _ambigresolver[iambig];
    }
  }

  void KalFit::makeTrack(TrkDef const& mytrk,TrkKalFit& myfit) {
    myfit._fit = TrkErrCode(TrkErrCode::fail);
// test if fitable
    if(fitable(mytrk)){
// create the hits. This also initializes T0
      makeHits(mytrk,myfit);
// Create the BaBar hit list, and fill it with these hits.  The BaBar list takes ownership
      std::vector<DetIntersection> detinter;
      TrkHotListFull* hotlist = new TrkHotListFull();
      for(std::vector<TrkStrawHit*>::iterator ihit=myfit._hits.begin();ihit!=myfit._hits.end();ihit++){
        TrkStrawHit* trkhit = *ihit;
	hotlist->append(trkhit);
// Also extract wall and gas intersection objects from each straw hit (active or not)
	DetIntersection wallinter;
	wallinter.delem = 0;
	wallinter.pathlen = trkhit->fltLen();
	DetIntersection gasinter;
	gasinter.delem = 0;
	gasinter.pathlen = trkhit->fltLen();
	if(mytrk.traj() != 0){
	  if(trkhit->wallElem().reIntersect(mytrk.traj(),wallinter))
	    detinter.push_back(wallinter);
	  if(trkhit->gasElem().reIntersect(mytrk.traj(),gasinter))
	    detinter.push_back(gasinter);
	} else {
	  if(trkhit->wallElem().reIntersect(&mytrk.helix(),wallinter))
	    detinter.push_back(wallinter);
	  if(trkhit->gasElem().reIntersect(&mytrk.helix(),gasinter))
	    detinter.push_back(gasinter);
	}
      }
// Create BaBar track
      myfit._trk = new TrkRecoTrk(_kalcon->defaultType(), 0, 0);
      assert(myfit._trk != 0);
      myfit._trk->setBField(bField());
      myfit._trk->resetT0(myfit._t0.t0(),myfit._t0.t0Err());
// create Kalman rep
      if(mytrk.traj() != 0)
	myfit._krep = new KalRep(mytrk.traj(), hotlist, detinter, myfit._trk, *_kalcon, PdtPid::electron);
      else
	myfit._krep = new KalRep(mytrk.helix(), hotlist, detinter, myfit._trk, *_kalcon, PdtPid::electron);
      assert(myfit._krep != 0);
      myfit._trk->setRep(myfit._krep);
// now fit
      fitTrack(myfit);
      if(_removefailed)myfit.removeFailed();
    }
  }

  void KalFit::addHits(TrkKalFit& myfit,const StrawHitCollection* straws, std::vector<hitIndex> indices) {
// there must be a valid Kalman fit to add hits to
    if(myfit._krep != 0 && myfit._fit.success()){
      ConditionsHandle<TrackerCalibrations> tcal("ignored");
      const Tracker& tracker = getTrackerOrThrow();
// setup to estimate initial flightlength
      Hep3Vector tdir;
      HepPoint tpos;
      double flt0 = findZFltlen(myfit,0.0);
      myfit._krep->referenceTraj()->getInfo(0.0,tpos,tdir);
      for(unsigned iind=0;iind<indices.size(); ++iind){
	size_t istraw = indices[iind]._index;
	const StrawHit& strawhit(straws->at(istraw));
	const Straw& straw = tracker.getStraw(strawhit.strawIndex());
	// estimate  initial flightlength
	double hflt = findZFltlen(myfit,straw.getMidPoint().z());
	// create the hit object.  Assume we're at the last iteration over added error
	TrkStrawHit* trkhit = new TrkStrawHit(strawhit,straw,istraw,myfit._t0,flt0,hflt,_herr.back(),_maxdriftpull);
	assert(trkhit != 0);
	if(indices[iind]._ambig != 0)trkhit->setAmbig(indices[iind]._ambig);
	// flag the added hit
	trkhit->setUsability(3);
// add the hit to the track and the fit
	myfit._krep->addHot(trkhit);
	myfit._hits.push_back(trkhit);
// create intersections for the material of this hit and add those to the track
	DetIntersection wallinter;
	if(trkhit->wallElem().reIntersect(myfit._krep->referenceTraj(),wallinter))
	  myfit._krep->addInter(wallinter);	
	DetIntersection gasinter;
	if(trkhit->gasElem().reIntersect(myfit._krep->referenceTraj(),gasinter))
	  myfit._krep->addInter(gasinter);	
      }
// refit the track
      myfit._krep->resetFit();
      fitTrack(myfit);
    }
  }

  void KalFit::fitTrack(TrkKalFit& myfit){
    // loop over external hit errors
    for(size_t iherr=0;iherr < _herr.size(); ++iherr){
// update the external hit errors.  This isn't strictly necessary on the 1st iteration.
      for(std::vector<TrkStrawHit*>::iterator itsh = myfit._hits.begin(); itsh != myfit._hits.end(); ++itsh){
	(*itsh)->setExtErr(_herr[iherr]);
      }
      // update t0, and propagate it to the hits
      double oldt0(myfit._t0.t0());
      myfit._nt0iter = 0;
      unsigned niter(0);
      bool changed(true);
//      myfit.fit();
      while(changed && niter < _kalcon->maxIterations()){
	changed = false;
	_ambigresolver[iherr]->resolveTrk(myfit);
	myfit._krep->resetFit();
	myfit.fit();
	if(! myfit._fit.success())break;
	if(_updatet0){
	  updateT0(myfit);
	  changed |= fabs(myfit._t0.t0()-oldt0) > _t0tol;
	  oldt0 = myfit._t0.t0();
	}
	// drop outlyers
	if(_weedhits){
	  myfit._nweediter = 0;
	  changed |= weedHits(myfit);
	}
	niter++;
      }
      if(myfit._krep != 0) myfit._krep->addHistory(myfit._fit,"KalFit");
    }
  }

  bool
  KalFit::fitable(TrkDef const& mytrk){
    return mytrk.strawHitIndices().size() >= _minnstraws;
  }

  void
  KalFit::initT0(TrkDef const& mytrk,TrkT0& t00) {
    // depending on the strategy, either compute T0 from the hits, or take it from the existing defintion directly
    if(_t0strategy == external){
      t00 = mytrk.trkT0();
    } else {
      // make an array of all the hit times, correcting for propagation delay
      const Tracker& tracker = getTrackerOrThrow();
      ConditionsHandle<TrackerCalibrations> tcal("ignored");
      unsigned nind = mytrk.strawHitIndices().size();
      std::vector<double> times;
      times.reserve(nind);
      // get flight distance of z=0 (for comparison)
      double t0flt = mytrk.helix().zFlight(0.0);
      // loop over strawhits
      for(unsigned iind=0;iind<nind;iind++){
	size_t istraw = mytrk.strawHitIndices()[iind]._index;
	const StrawHit& strawhit(mytrk.strawHitCollection()->at(istraw));
	const Straw& straw = tracker.getStraw(strawhit.strawIndex());
	// assume the average drift time is half the maximum drift distance.  This is a poor approximation, but good enough for now
	D2T d2t;
	static CLHEP::Hep3Vector zdir(0.0,0.0,1.0);
	tcal->DistanceToTime(straw.index(),0.5*straw.getRadius(),zdir,d2t);
	// compute initial flightlength from helix and hit Z
	double hflt = mytrk.helix().zFlight(straw.getMidPoint().z()) - t0flt;
	// estimate the time the track reaches this hit from z=0, assuming speed-of-light travel along the helix. Should use actual
	// speed based on momentum and assumed particle species FIXME!!!
	double tprop = hflt/_vlight;
	// estimate signal propagation time on the wire assuming the middle (average)
	double vwire = tcal->SignalVelocity(straw.index());
	double teprop = straw.getHalfLength()/vwire;
	// correct the measured time for these effects: this gives the aveage time the particle passed this straw, WRT
	// when the track crossed Z=0
	double htime = strawhit.time() - tprop - teprop - d2t._tdrift;
	times.push_back(htime);
      }
      // different strategies for t0 estimate
      if(_t0strategy == median){
	// find the median hit time
	std::sort(times.begin(),times.end());
	unsigned imed = times.size()/2;
	// deal with even/odd # of hits separately
	if(times.size() == 2*imed)
	  t00._t0 = 0.5*(times[imed-1]+times[imed]);
	else
	  t00._t0 = times[imed];
	// estimate the error using the range
	double tmax = (times.back()-times.front());
	t00._t0err = _t0errfac*tmax/sqrt(12*nind);
      } else if(_t0strategy == histogram) {
      }
    }
  }
  
  void
  KalFit::makeHits(TrkDef const& mytrk, TrkKalFit& myfit) {
// first, find t0.  This can come from the hits or from the track definition
    TrkT0 t00;
    initT0(mytrk,t00);
    myfit._t0 = t00;
    myfit._t00 = t00;
    const Tracker& tracker = getTrackerOrThrow();
    unsigned nind = mytrk.strawHitIndices().size();
    double flt0 = mytrk.helix().zFlight(0.0);
    for(unsigned iind=0;iind<nind;iind++){
      size_t istraw = mytrk.strawHitIndices()[iind]._index;
      const StrawHit& strawhit(mytrk.strawHitCollection()->at(istraw));
      const Straw& straw = tracker.getStraw(strawhit.strawIndex());
      double fltlen = mytrk.helix().zFlight(straw.getMidPoint().z());
    // create the hit object.  Start with the 1st additional error for anealing
      TrkStrawHit* trkhit = new TrkStrawHit(strawhit,straw,istraw,t00,flt0,fltlen,_herr.front(),_maxdriftpull);
      assert(trkhit != 0);
    // set the initial ambiguity based on the input
      trkhit->setAmbig(mytrk.strawHitIndices()[iind]._ambig);
    // refine the flightlength, as otherwise hits in the same plane are at exactly the same flt, which can cause problems
      const TrkDifTraj* dtraj = &mytrk.helix();
      if(mytrk.traj() != 0)dtraj = mytrk.traj();
      TrkErrCode pstat = trkhit->updatePoca(dtraj);
      if(pstat.failure()){
        trkhit->setActivity(false);
      }
      myfit._hits.push_back(trkhit);
    }
 // sort the hits by flightlength
    std::sort(myfit._hits.begin(),myfit._hits.end(),fltlencomp());
  }
  
  bool
  KalFit::updateT0(TrkKalFit& myfit){
    bool retval(false);
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
// need to have a valid fit
    if(myfit._krep->fitValid()){
// find the global fltlen associated with z=0.  This should be a piectraj function, FIXME!!!
      double flt0 = findZFltlen(myfit,0.0);
// find hits
      std::vector<double> hitst0; // store t0, to allow outlyer removal
      double t0sum(0.0);
      double t0sum2(0.0);
      for(std::vector<TrkStrawHit*>::iterator ihit= myfit._hits.begin();ihit != myfit._hits.end(); ihit++){
        TrkStrawHit* hit = *ihit;
        if(hit->isActive() && hit->poca()!= 0 && hit->poca()->status().success()){
// copy the seed
          static TrkSimpTraj* straj = myfit._krep->seed()->clone();
// find the hit site in the rep
          const KalHit* hitsite = myfit._krep->findHotSite(hit);
// set helix to the local parameters EXCLUDING THIS HIT
          if(hitsite != 0 && myfit._krep->smoothedTraj(hitsite,straj)){
            TrkPoca poca(*straj,hit->fltLen(),*(hit->hitTraj()),hit->hitLen());
            if(poca.status().success()){
// Sign doca by the ambiguity
	      double doca = poca.doca()*hit->ambig();
// restrict the range, symmetrically to avoid bias
	      double rad = hit->straw().getRadius();
              if(doca > _mint0doca && doca < rad-_mint0doca){
// translate the DOCA into a time
		D2T d2t;
		tcal->DistanceToTime(hit->straw().index(),doca,straj->direction(poca.flt1()),d2t);
// particle transit time to this hit from z=0.  This assumes beta=1, FIXME!!!
                double tflt = (hit->fltLen()-flt0)/_vlight;
// electronic signal transit time from the POCA to the primary straw end.
// This calculation corresponds tot he conventions that the time is measured WRT the far end of the wire as defined
// by the wire direction, and that the hit trajectory starts in the middle of the straw, going in the wire direction.
		double teprop = (hit->straw().getHalfLength()-poca.flt2())/tcal->SignalVelocity(hit->straw().index()); 
// t0 = Time difference between the hit time and the sum of all propagation
                double hitt0 = hit->time() - d2t._tdrift - tflt - teprop;
                hitst0.push_back(hitt0);
		t0sum += hitt0;
		t0sum2 += hitt0*hitt0;
              }
            }
          }
        }
      }
      if(hitst0.size() >1){
// find the median, then average hits in a window around that.
	std::sort(hitst0.begin(),hitst0.end());
	unsigned imed = hitst0.size()/2;
	double t0med;
	if(hitst0.size() == 2*imed)
	// even
	  t0med = 0.5*(hitst0[imed-1]+hitst0[imed]);
	else
	// odd
	  t0med = hitst0[imed];
// initial sigma from all hits
	double t0mean = t0sum/hitst0.size();
	double t02 = t0sum2/hitst0.size();
	double t0sig = sqrt(max(t02 - t0mean*t0mean,0.0));
        double t0(t0med);
        double t0err(-1.0);
	bool changed(true);
// iterate until the set of used hits doesn't change
	std::vector<bool> used(hitst0.size(),true);
	unsigned niter(0);
	unsigned nactive(hitst0.size());
	while(changed && niter < 10 && nactive > 1){
	  niter++;
	  changed = false;
	  t0sum = t0sum2 = 0.0;
	  nactive = 0;
	  for(unsigned ihit=0;ihit<hitst0.size();ihit++){
	    bool useit = fabs(hitst0[ihit]-t0) < _t0nsig*t0sig;
	    changed |= useit != used[ihit];
	    used[ihit] = useit;
	    if(useit){
	      t0sum += hitst0[ihit];
	      t0sum2 += hitst0[ihit]*hitst0[ihit];
	      ++nactive;
	    }
	  }
	  if(nactive > 1){
	    t0 = t0sum/nactive;
	    double t02 = t0sum2/nactive;
	    double t0sig = sqrt(max(t02 - t0*t0,0.0));
	    t0err = _t0errfac*t0sig/sqrt(nactive-1);
	  }
	}
// reset t0
        myfit._t0.setT0(t0,t0err);
	if(myfit._trk != 0)myfit._trk->resetT0(t0,t0err);
// reset all the hit times
        for(std::vector<TrkStrawHit*>::iterator ihit= myfit._hits.begin();ihit != myfit._hits.end(); ihit++){
          (*ihit)->updateT0(myfit._t0,flt0);
        }
        retval = true;
      }
    }
    return retval;
  }

  bool
  KalFit::weedHits(TrkKalFit& myfit) {
    // Loop over HoTs and find HoT with largest contribution to chi2.  If this value
    // is greater than some cut value, deactivate that HoT and reFit
    bool retval(false);
    double worst = -1.;
    TrkHitOnTrk* worstHot = 0;
    TrkHotList* hots = myfit._krep->hotList();
    for (TrkHotList::nc_hot_iterator iHot = hots->begin(); iHot != hots->end(); ++iHot) {
      if (iHot->isActive()) {
        double resid, residErr;
        if(iHot->resid(resid, residErr, true)){
          double value = fabs(resid/residErr);
          if (value > _maxhitchi && value > worst) {
            worst = value;
            worstHot = iHot.get();
          }
        }
      }
    }
    if(0 != worstHot){
      retval = true;
      worstHot->setActivity(false);
      worstHot->setUsability(-5);
      myfit.fit();
      myfit._krep->addHistory(myfit._fit, "HitWeed");
      // Recursively iterate
      myfit._nweediter++;
      if (myfit._fit.success() && myfit._nweediter < _maxweed ) {
        retval |= weedHits(myfit);
      }
    }
    return retval;
  }

  const BField*
  KalFit::bField() {
    if(_bfield == 0){
// create a wrapper around the mu2e nominal DS field
      GeomHandle<BFieldConfig> bfconf;
      _bfield=new BFieldFixed(bfconf->getDSUniformValue());
      assert(_bfield != 0);
    }
    return _bfield;
  }

  double KalFit::findZFltlen(const TrkKalFit& myfit,double zval) {
    double loclen;
    double zflt = zval/myfit._krep->traj().direction(0.0).z();
    double dz(10.0);
    unsigned niter = 0;
    while(fabs(dz) > 1.0 && niter < _kalcon->maxIterations() ) {
      const HelixTraj* helix = dynamic_cast<const HelixTraj*>(myfit._krep->localTrajectory(zflt,loclen));
      zflt += helix->zFlight(0.0)-loclen;
      dz = myfit._krep->traj().position(zflt).z();
      niter++;
    }
    return zflt;
  }
}


