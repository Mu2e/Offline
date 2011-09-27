//
// Class to perform BaBar Kalman fit
//
// $Id: KalFit.cc,v 1.14 2011/09/27 21:49:09 mu2ecvs Exp $
// $Author: mu2ecvs $ 
// $Date: 2011/09/27 21:49:09 $
//

// the following has to come before other BaBar includes
#include "BaBar/BaBar.hh"
#include "KalmanTests/inc/KalFit.hh"
//geometry
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
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
// drift velocity should come from a service FIXME!!!  
  const double KalFit::_vdrift = 0.05; // 50 um/nsec
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
    _material(pset.get<bool>("material",false)),
    _ambigflip(pset.get<bool>("ambigflip",false)),
    _weedhits(pset.get<bool>("weedhits",true)),
    _updatet0(pset.get<bool>("updateT0",true)),
    _removefailed(pset.get<bool>("RemoveFailedFits",true)),
    _t0tol(pset.get<double>("t0Tolerance",1.0)),
    _maxhitchi(pset.get<double>("maxhitchi",5.0)),
    _maxiter(pset.get<unsigned>("maxiter",3)),
    _mingap(pset.get<double>("mingap",0.1)),
    _minnstraws(pset.get<unsigned>("minnstraws",20)),
    _minndof(pset.get<unsigned>("minNDOF",15)),
    _maxweed(pset.get<unsigned>("maxweed",10)),
    _herr(pset.get<double>("hiterr",0.1)),
    _ssmear(pset.get<double>("seedsmear",1e6)),
    _t0errfac(pset.get<double>("t0ErrorFactor",1.2))
  {
      _kalcon = new KalContext;
      _kalcon->setBendSites(_fieldcorr);
      _kalcon->setMaterialSites(_material);
      _kalcon->setForbidAmbigFlips(_ambigflip); //false: free left-rigth ambiguity, true: will be fixed from sim
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
      _kalcon->setDefaultType(PdtPid::electron); // by default, fit electrons
  }

  KalFit::~KalFit(){
    delete _kalcon;
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
	DetIntersection gasinter;
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

  void KalFit::addHits(TrkKalFit& myfit,const StrawHitCollection* straws, std::vector<size_t> indices) {
// there must be a valid Kalman fit to add hits to
    if(myfit._krep != 0 && myfit._fit.success()){
      ConditionsHandle<TrackerCalibrations> tcal("ignored");
      const Tracker& tracker = getTrackerOrThrow();
// setup to estimate initial flightlength
      Hep3Vector tdir;
      HepPoint tpos;
      myfit._krep->referenceTraj()->getInfo(0.0,tpos,tdir);
      for(unsigned iind=0;iind<indices.size(); ++iind){
	size_t istraw = indices[iind];
	const StrawHit& strawhit(straws->at(istraw));
	const Straw& straw = tracker.getStraw(strawhit.strawIndex());
	// estimate  initial flightlength
	double hflt = straw.getMidPoint().z()/tdir.z();
	// estimate the time the track reaches this hit, assuming speed-of-light travel along the helix. Should use actual
	// speed based on momentum and assumed particle species FIXME!!!
	double tprop = hflt/_vlight;
	double hitt0 = myfit._t0.t0() + tprop;
	// create the hit object
	TrkStrawHit* trkhit = new TrkStrawHit(strawhit,straw,istraw,hitt0,myfit._t0.t0Err(),_herr);
	// flag the added hit
	trkhit->setUsability(3);
	assert(trkhit != 0);
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
    // fit the track
    myfit.fit();
    // update t0, and propagate it to the hits
    double oldt0(-1e8);
    myfit._nt0iter = 0;
    while(_updatet0 && myfit._fit.success() && fabs(myfit._t0.t0()-oldt0) > _t0tol && 
	myfit._nt0iter < _kalcon->maxIterations()){
      oldt0 = myfit._t0.t0();
      if(updateT0(myfit)){
	myfit._krep->resetFit();
	myfit.fit();
	myfit._nt0iter++;
      } else
	break;
      // drop outlyers
      if(_weedhits){
	myfit._nweediter = 0;
	weedHits(myfit);
      }
    }
  }

  
  bool
  KalFit::fitable(TrkDef const& mytrk){
    return mytrk.strawHitIndices().size() >= _minnstraws;
  }
  
  void
  KalFit::makeHits(TrkDef const& mytrk, TrkKalFit& myfit) {
    const Tracker& tracker = getTrackerOrThrow();
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
    unsigned nind = mytrk.strawHitIndices().size();
    std::vector<double> times;
    times.reserve(nind);
    double flt0 = mytrk.helix().zFlight(0.0);
    for(unsigned iind=0;iind<nind;iind++){
      unsigned istraw = mytrk.strawHitIndices()[iind];
      const StrawHit& strawhit(mytrk.strawHitCollection()->at(istraw));
      const Straw& straw = tracker.getStraw(strawhit.strawIndex());
    // compute initial flightlength from helix and hit Z
      double hflt = mytrk.helix().zFlight(straw.getMidPoint().z());
    // estimate the time the track reaches this hit from z=0, assuming speed-of-light travel along the helix. Should use actual
    // speed based on momentum and assumed particle species FIXME!!!
      double tprop = hflt/_vlight;
      double hitt0 = mytrk.trkT0().t0() + tprop;
    // subtract the propagation time and the average wire signal delay when computing hit time
      double vwire = tcal->SignalVelocity(straw.index());
      double htime = strawhit.time() - tprop - straw.getHalfLength()/vwire;
      times.push_back(htime);
    // create the hit object
      TrkStrawHit* trkhit = new TrkStrawHit(strawhit,straw,istraw,hitt0,mytrk.trkT0().t0Err(),_herr);
      assert(trkhit != 0);
    // refine the flightlength, as otherwise hits in the same plane are at exactly the same flt, which can cause problems
      const TrkDifTraj* dtraj = &mytrk.helix();
      if(mytrk.traj() != 0)dtraj = mytrk.traj();
      TrkErrCode pstat = trkhit->updatePoca(dtraj);
      if(pstat.failure()){
        trkhit->setActivity(false);
      }
      myfit._hits.push_back(trkhit);
    }
  // if the initial t0error was < 0, compute t0 and override
    if(mytrk.trkT0().t0Err() < 0.0 && nind > 0){
  // find the median hit time
      std::sort(times.begin(),times.end());
      unsigned imed = times.size()/2;
      double tmed;
      if(times.size() == 2*imed)
	tmed = 0.5*(times[imed-1]+times[imed]);
      else
	tmed = times[imed];
  // assuming a flat drift time means correcting by half the maximum drift time
      double tmax = myfit._hits[0]->straw().getRadius()/_vdrift;
      double t0 = tmed - 0.5*tmax;
  // estimate the error using the same assumption
      double t0err = _t0errfac*tmax/sqrt(12*nind);
  // save this as the initial T0 value (and the current t0 value)
      myfit.setT00(TrkT0(t0,t0err));
      myfit.setT0(TrkT0(t0,t0err));
    } else {
  // initialize t0 directly from the tracking object
      myfit.setT00(mytrk.trkT0());
      myfit.setT0(mytrk.trkT0());
    }
  // update the hits
    for(std::vector<TrkStrawHit*>::iterator ihit= myfit._hits.begin();ihit != myfit._hits.end(); ihit++){
      double hitt0 = myfit._t0.t0() + ((*ihit)->fltLen() -flt0)/_vlight;
      (*ihit)->updateT0(hitt0,myfit._t0.t0Err());
    }
  // sort the hits by flightlength
    std::sort(myfit._hits.begin(),myfit._hits.end(),fltlencomp());
  }
  
  bool
  KalFit::updateT0(TrkKalFit& myfit){
    bool retval(false);
// need to have a valid fit
    if(myfit._krep->fitValid()){
// find the global fltlen associated with z=0.  This should be a piectraj function, FIXME!!!
      double loclen;
      double flt0 = 0.0;
      double dz(10.0);
      unsigned niter = 0;
      while(fabs(dz) > 1.0 && niter < _kalcon->maxIterations() ) {
        const HelixTraj* helix = dynamic_cast<const HelixTraj*>(myfit._krep->localTrajectory(flt0,loclen));
        flt0 += helix->zFlight(0.0)-loclen;
        dz = myfit._krep->traj().position(flt0).z();
        niter++;
      }
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
// DOCA can't be more than the straw radius!
              double doca = std::min(fabs(poca.doca()),hit->straw().getRadius());
// require a minimum doca to avoid ambiguity bias.  mindoca shoudl be a parameter, FIXME!!!
              static double mindoca(0.0);
              if(doca > mindoca){
// propagation time to this hit from z=0.  This assumes beta=1, FIXME!!!
                double tflt = (hit->fltLen()-flt0)/_vlight;
// drift time of this hit (plus t0)
                double tdrift = hit->time() - tflt;
// t0 = Time difference between the drift time and the DOCA time.  sign of DOCA is irrelevant here.
                double hitt0 = tdrift - doca/_vdrift;
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
// TrkT0 window should be a parameter, FIXME!!!!
        double nsig(2.5);
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
	    bool useit = fabs(hitst0[ihit]-t0) < nsig*t0sig;
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
// reset all the hit times
        for(std::vector<TrkStrawHit*>::iterator ihit= myfit._hits.begin();ihit != myfit._hits.end(); ihit++){
          TrkStrawHit* hit = *ihit;
// correct for flightlength.  Again assumes beta=1, FIXME!!!
          double hitt0 = t0 + (hit->fltLen()-flt0)/_vlight;
          hit->updateT0(hitt0,t0err);
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
      GeomHandle<BFieldManager> bfMgr;
      _bfield=new BFieldFixed(bfMgr->getDSUniformValue());
      assert(_bfield != 0);
    }
    return _bfield;
  }
}

