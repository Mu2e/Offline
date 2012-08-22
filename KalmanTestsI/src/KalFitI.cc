//
// Class to perform BaBar Kalman fit
//
// $Id: KalFitI.cc,v 1.1 2012/08/22 17:30:37 tassiell Exp $
// $Author: tassiell $ 
// $Date: 2012/08/22 17:30:37 $
//

// the following has to come before other BaBar includes
#include "BaBar/BaBar.hh"

#include "KalmanTestsI/inc/KalFitI.hh"

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
#include "DchGeom/DchDetector.hh"
#include "DetectorModel/DetSet.hh"
#include "ITrackerGeom/inc/ITracker.hh"

#include "TMath.h"
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <math.h>


using namespace std; 

namespace mu2e 
{
// convert speed of light in mm/nsec
  const double KalFitI::_vlight = CLHEP::c_light;
// drift velocity should come from a service FIXME!!!  
  const double KalFitI::_vdrift = 0.035; // 50 um/nsec
// comparison functor for ordering hits
  /*struct fltlencomp : public binary_function<TrkStrawHit*, TrkStrawHit*, bool> {
          bool operator()(TrkStrawHit* x, TrkStrawHit* y) { return x->fltLen() < y->fltLen(); }
  };*/
  struct fltlencomp : public binary_function<TrkStrawHit*, TrkStrawHit*, bool> {
    fltlencomp(TrkFitDirection::FitDirection fdir=TrkFitDirection::downstream) : _fdir(fdir) {}
    //fltlencomp(KalFit::fitDirection fdir=KalFit::downstream) : _fdir(fdir) {}
    bool operator()(TrkStrawHit* x, TrkStrawHit* y) { 
      return _fdir == TrkFitDirection::downstream ? x->fltLen() < y->fltLen() : y->fltLen() < x->fltLen() ;
      //return _fdir == KalFit::downstream ? x->fltLen() < y->fltLen() : y->fltLen() < x->fltLen() ;
    }
    TrkFitDirection::FitDirection _fdir;
    //KalFit::fitDirection _fdir;
  };


// construct from a parameter set  
  KalFitI::KalFitI(fhicl::ParameterSet const& pset) : _bfield(0),
    _debug(pset.get<int>("debugLevel",0)),
    _fieldcorr(pset.get<bool>("fieldCorrections",false)),
    _material(pset.get<bool>("material",true)),
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
    //_herr(pset.get<double>("hiterr",0.2)),
    _herr(pset.get< vector<double> >("hiterr")),
    _ssmear(pset.get<double>("seedsmear",1e6)),
    _t0errfac(pset.get<double>("t0ErrorFactor",1.2)),
    _mint0doca(pset.get<double>("minT0DOCA",-0.2)),
    _maxdriftpull(pset.get<double>("maxDriftPull",10)),
    _t0nsig(pset.get<double>("t0window",2.5)),
    _tpart((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
    _t0strategy((t0Strategy)pset.get<int>("t0strategy",median)),
    _exacthitturn(pset.get<bool>("exacthitturn",false)),
    _ambigstrategy(pset.get< vector<int> >("ambiguityStrategy"))
  {
      // make sure we have at least one entry for additional errors
      if(_herr.size() <= 0) throw cet::exception("RECO")<<"mu2e::KalFitI: no hit errors specified" << endl;
      if(_herr.size() != _ambigstrategy.size()) throw cet::exception("RECO")<<"mu2e::KalFitI: inconsistent ambiguity resolution" << endl;
      _kalcon = new KalContext;
      _kalcon->setBendSites(_fieldcorr);
      _kalcon->setMaterialSites(_material);
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
      // _kalcon->setDefaultType(_tpart);
      // construct the ambiguity resolvers
            for(size_t iambig=0;iambig<_ambigstrategy.size();++iambig){
              switch (_ambigstrategy[iambig] ){
                case KalFit::fixedambig: default:
                  _ambigresolver.push_back(new FixedAmbigResolver(pset));
                  break;
                case KalFit::pocaambig:
                  _ambigresolver.push_back(new PocaAmbigResolver(pset));
                  break;
                case KalFit::hitambig:
                  _ambigresolver.push_back(new HitAmbigResolver(pset));
                  break;
                case KalFit::panelambig:
                  _ambigresolver.push_back(new PanelAmbigResolver(pset));
                  break;
              }
            }

    }

  KalFitI::~KalFitI(){
    delete _kalcon;
    for(size_t iambig=0;iambig<_ambigresolver.size();++iambig){
      delete _ambigresolver[iambig];
    }
  }

  void KalFitI::makeTrack(TrkDef& mytrk,TrkKalFit& myfit) {
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
      }

      HelixTraj*  seed=mytrk.helix().clone();
      //find material intersection before fitting
      double helixturn=2*TMath::Pi()/seed->omega()/seed->cosDip();
      hotlist->sort();
      double range[2];
      range[0] = hotlist->startFoundRange();
      range[0] = (range[0]>0?0.:-helixturn);
      range[1] = hotlist->endFoundRange()+0.001;
      seed->setFlightRange(range);
      //seed.printAll(std::cout);
      TrkDifPieceTraj* reftraj = new TrkDifPieceTraj(seed,range[0],range[1]);
      //reftraj->printAll(std::cout);
      const DetSet* trkmodel= &(DchDetector::GetInstance()->dchSet());
      trkmodel->intersection(detinter,reftraj,range);
      //reftraj->printAll(std::cout);
      delete reftraj;
      if(_debug>1) std::cout<<"created detector intersect for range "<<range[0]<<" "<<range[1]<<" size="<<detinter.size()<<std::endl;

      
// create Kalman rep
      if(mytrk.traj() != 0)
        myfit._krep = new KalRep(mytrk.traj(), hotlist, detinter,  *_kalcon, mytrk.particle());
      else
	myfit._krep = new KalRep(mytrk.helix(), hotlist, detinter, *_kalcon, mytrk.particle());
      assert(myfit._krep != 0);
      myfit._krep->setT0(myfit._t0);

      if(_debug>5) {
	std::cout<<"print kalrep"<<std::endl;
	myfit._krep->printAll(std::cout);
      }
// now fit
      fitTrack(myfit);
      //      if(_removefailed)myfit.removeFailed();
    }
  }

  void KalFitI::addHits(TrkDef& mytrk,TrkKalFit& myfit, std::vector<hitIndex> indices,bool active) {
// there must be a valid Kalman fit to add hits to
    if(myfit._krep != 0 && myfit._fit.success()){
      mytrk.setTraj(&myfit._krep->pieceTraj());
      const Tracker& tracker = getTrackerOrThrow();
      const mu2e::ITracker &itracker = static_cast<const mu2e::ITracker&>( tracker );    
      CellGeometryHandle *itwp = itracker.getCellGeometryHandle();
// setup to estimate initial flightlength
      Hep3Vector tdir;
      HepPoint tpos;
      double flt0 = findZFltlen(myfit,0.0);
      myfit._krep->referenceTraj()->getInfo(0.0,tpos,tdir);
      for(unsigned iind=0;iind<indices.size(); ++iind){
	size_t istraw = indices[iind]._index;
	const StrawHit& strawhit(mytrk.strawHitCollection()->at(istraw));
	itwp->SelectCellDet(strawhit.strawIndex().asUint());
	double zmid=itwp->GetCellCenter().z();
	const Straw& straw = tracker.getStraw(strawhit.strawIndex());
	// estimate  initial flightlength
	double hflt = _exacthitturn?(flt0+_hitflt[iind]):findZFltlen(myfit,zmid);
	// create the hit object
//	TrkStrawHit* trkhit = new TrkStrawHit(strawhit,straw,istraw,myfit._t0,flt0,hflt,_herr,_maxdriftpull);
        TrkStrawHit* trkhit = new TrkStrawHit(strawhit,straw,istraw,myfit._t0,flt0,hflt,_herr.back(),_maxdriftpull);
	assert(trkhit != 0);
	if(indices[iind]._ambig != 0)trkhit->setAmbig(indices[iind]._ambig);
	fixHitTurn(mytrk,trkhit);
	// flag the added hit
	trkhit->setUsability(3);
	trkhit->setActivity(active);
// add the hit to the track and the fit
	myfit._krep->addHot(trkhit);
	myfit._hits.push_back(trkhit);
      }
// refit the track
      myfit._krep->resetFit();
      fitTrack(myfit);
    }
  }

  void KalFitI::addHitsUnique(TrkDef& mytrk,TrkKalFit& myfit, std::vector<hitIndex> indices,bool active) {
    std::vector<hitIndex> missed;
    for(unsigned iind=0;iind<indices.size(); ++iind){
      size_t istraw = indices[iind]._index;
      const StrawHit& strawhit(mytrk.strawHitCollection()->at(istraw));
      std::vector<TrkStrawHit*>::iterator ifnd = find_if(myfit._hits.begin(),myfit._hits.end(),FindTrkStrawHit(strawhit));
      if(ifnd == myfit._hits.end()){
	missed.push_back(istraw);
      }
    }
    addHits(mytrk,myfit,missed,active);
  }

  void KalFitI::fitTrack(TrkKalFit& myfit){
    // loop over external hit errors
    for(size_t iherr=0;iherr < _herr.size(); ++iherr){
      // update the external hit errors.  This isn't strictly necessary on the 1st iteration.
      for(std::vector<TrkStrawHit*>::iterator itsh = myfit._hits.begin(); itsh != myfit._hits.end(); ++itsh){
        (*itsh)->setExtErr(_herr[iherr]);
      }


    // fit the track
    myfit._fit = TrkErrCode::succeed;
    _ambigresolver[iherr]->resolveTrk(myfit);
    myfit.fit();
    std::cout<<"chi2 "<<myfit._krep->chisquared(trkIn)<<" t0 "<<myfit._t0.t0()<<" ndof "<<myfit._krep->nDof()<<std::endl;
    // update t0, and propagate it to the hits
    double oldt0(-1e8);
    myfit._nt0iter = 0;
    unsigned niter(0);
    bool changed=true;
    while(changed && myfit._fit.success() && niter < _kalcon->maxIterations()){
      changed = false;
      if(_updatet0 && myfit._nt0iter < _kalcon->maxIterations() && updateT0(myfit) && fabs(myfit._t0.t0()-oldt0) > _t0tol  ){
	oldt0 = myfit._t0.t0();
        _ambigresolver[iherr]->resolveTrk(myfit);
	myfit._krep->resetFit();
	myfit.fit();
	std::cout<<"chi2 "<<myfit._krep->chisquared(trkIn)<<" t0 "<<myfit._t0.t0()<<" ndof "<<myfit._krep->nDof()<<std::endl;
	myfit._nt0iter++;
	changed = true;
      }
      // drop outlyers
      if(_weedhits){
	myfit._nweediter = 0;
	changed |= weedHits(myfit);
	std::cout<<"afterweedchi2 "<<myfit._krep->chisquared(trkIn)<<" t0 "<<myfit._t0.t0()<<" ndof "<<myfit._krep->nDof()<<std::endl;
      }
      niter++;
    }
    if(myfit._krep != 0) myfit._krep->addHistory(myfit._fit,"KalFitI");

    }//close loop on _herr

  }

  
  bool
  KalFitI::fitable(TrkDef const& mytrk){
    return mytrk.strawHitIndices().size() >= _minnstraws;
  }
  
  void
  KalFitI::initT0(TrkDef const& mytrk,TrkT0& t00) {
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
	// assume a constant drift velocity means the average drift time is half the maximum drift time
	double tdrift = 0.5*straw.getRadius()/_vdrift;
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
	double htime = strawhit.time() - tprop - teprop - tdrift;
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
  KalFitI::makeHits(TrkDef const& mytrk, TrkKalFit& myfit) {
// first, find t0.  This can come from the hits or from the track definition
    TrkT0 t00;
    initT0(mytrk,t00);
    const Tracker& tracker = getTrackerOrThrow();
    const mu2e::ITracker &itracker = static_cast<const mu2e::ITracker&>( tracker );    
    CellGeometryHandle *itwp = itracker.getCellGeometryHandle();

    unsigned nind = mytrk.strawHitIndices().size();
    double flt0 = mytrk.helix().zFlight(0.0);
    if(_t0strategy == external){
      t00._t0+=(flt0-_flt0)/CLHEP::c_light;
    }
    myfit._t0 = t00;
    myfit._t00 = t00;
    
    for(unsigned iind=0;iind<nind;iind++){
      size_t istraw = mytrk.strawHitIndices()[iind]._index;
      const StrawHit& strawhit(mytrk.strawHitCollection()->at(istraw));
      itwp->SelectCellDet(strawhit.strawIndex().asUint());
      double zmid=itwp->GetCellCenter().z();
      const Straw& straw = tracker.getStraw(strawhit.strawIndex());

      //estimate track length
      double fltlen = _exacthitturn?(_flt0+_hitflt[iind]):mytrk.helix().zFlight(zmid);
    // create the hit object
//      TrkStrawHit* trkhit = new TrkStrawHit(strawhit,straw,istraw,t00,flt0,fltlen,_herr,_maxdriftpull);
      TrkStrawHit* trkhit = new TrkStrawHit(strawhit,straw,istraw,t00,flt0,fltlen,_herr.front(),_maxdriftpull);
      assert(trkhit != 0);
      if(mytrk.strawHitIndices()[iind]._ambig != 0)trkhit->setAmbig(mytrk.strawHitIndices()[iind]._ambig);
    // refine the flightlength, as otherwise hits in the same plane are at exactly the same flt, which can cause problems
      fixHitTurn(mytrk,trkhit);
      myfit._hits.push_back(trkhit);
    }
 // sort the hits by flightlength
    std::sort(myfit._hits.begin(),myfit._hits.end(),fltlencomp(mytrk.fitdir().fitDirection()));
    //std::sort(myfit._hits.begin(),myfit._hits.end(),fltlencomp());
  }
  
  bool KalFitI::fixHitTurn(const TrkDef& mytrk,TrkStrawHit* trkhit){
    const TrkDifTraj* dtraj = &mytrk.helix();
    if(mytrk.traj() != 0)dtraj = mytrk.traj();
    
    //    const TrkDifTraj* dtraj = (mytrk.traj() != 0)?mytrk.traj():(&mytrk.helix());
    TrkErrCode pstat = trkhit->updatePoca(dtraj);
    //    double doca1=trkhit->poca()->doca();
    //    std::cout<<"t0 "<<trkhit->hitT0()<<" time "<<trkhit->time()<<" r "<<trkhit->driftRadius()
    //	     <<" doca "<<trkhit->poca()->doca()<<" res "<<trkhit->residual()<<std::endl;
    if(!_exacthitturn&&!pstat.failure()){
      const Tracker& tracker = getTrackerOrThrow();
      const mu2e::ITracker &itracker = static_cast<const mu2e::ITracker&>( tracker );    
      CellGeometryHandle *itwp = itracker.getCellGeometryHandle();
      itwp->SelectCellDet(trkhit->strawHit().strawIndex().asUint());
      double zmid=itwp->GetCellCenter().z();
      double dz=fabs(itwp->GetCellDirection().z()*itwp->GetCellHalfLength());

      const HelixTraj& seed=mytrk.helix();
      double helixzturn=fabs(2*TMath::Pi()/seed.omega()*seed.tanDip());
      double helixturn=helixzturn/seed.sinDip();
      double doca=fabs(trkhit->poca()->doca());
      double fltlen=trkhit->poca()->flt1();
      double hitlen=trkhit->poca()->flt2();
      double zstraw=hitlen;
      double safetyfactor=1.1;
      // int nturn1=int(fabs((zstraw+trkhit->straw().getHalfLength())/helixzturn)*safetyfactor);
      // int nturn2=int(fabs((trkhit->straw().getHalfLength()-zstraw)/helixzturn)*safetyfactor);
      int nturn1=ceil((zmid-dz*safetyfactor*safetyfactor-zstraw)/helixzturn);
      int nturn2=floor((zmid+dz*safetyfactor*safetyfactor-zstraw)/helixzturn);
      if(nturn1*nturn2>0) doca=1e12;
      int nmin=0;
      if(_debug>1) 
	cout<<trkhit->index()<<" can be "<<nturn1<<" "<<nturn2<<" from doca "<<doca
	    <<" fltlen "<<fltlen<<" hitlen "<<hitlen<<" helixzturn "<<helixzturn<<" z+-dzhalf wire legth "<<zmid<<" "<<dz<<endl;
      for(int i=nturn1;i<=nturn2;i++){
	if(i==0) continue;
	trkhit->setHitLen(hitlen+i*helixzturn);
	trkhit->setFltLen(fltlen+i*helixturn);
	TrkErrCode pstat = trkhit->updatePoca(dtraj);
	if(!pstat.failure()){
	  double doca1=fabs(trkhit->poca()->doca());
	  if(doca1<doca&&
	     fabs(trkhit->poca()->flt2()-zmid)<dz*safetyfactor){nmin=i;doca=doca1;} 
	  if(_debug>1) 
	    cout<<i<<" doca "<<doca1<<" fltlen "<<fltlen+i*helixturn
		<<" hitlen "<<hitlen+i*helixzturn
		<<" obtained "<<trkhit->poca()->flt1()<<" "<<trkhit->poca()->flt2()<<" min "<<nmin<<endl;
	}
      }
      //set to minimum value(can be taken from cache)
      if(nturn2-nturn1>=0||nmin!=0){
	trkhit->setHitLen(hitlen+nmin*helixzturn);
	trkhit->setFltLen(fltlen+nmin*helixturn);
	pstat = trkhit->updatePoca(dtraj); 
      }
      if(nmin!=0&&!pstat.failure()){
	trkhit->setActivity(true);
	trkhit->setUsability(3);
	return true;
      }
      //if(nmin!=0) cout<<"was selected wrong turn"<<endl;
    }
    if(pstat.failure()){
      trkhit->setActivity(false);
    }
    return false;
  }
  
  bool
  KalFitI::updateT0(TrkKalFit& myfit){
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
	      double rad = hit->straw().getRadius();
// sign doca by the ambiguity.  Restrict to the physical maximum
              double doca = std::min(poca.doca()*hit->ambig(),rad);
// restrict the range, symmetrically to avoid bias
              if(doca > _mint0doca){
// propagation time to this hit from z=0.  This assumes beta=1, FIXME!!!
                double tflt = (hit->fltLen()-flt0)/_vlight;
		double vwire = tcal->SignalVelocity(hit->straw().index());
		double eprop = (hit->straw().getHalfLength()-poca.flt2())/vwire; 
		// drift time of this hit (plus t0)
		double tdrift = hit->time() - tflt - eprop;
// t0 = Time difference between the drift time and the DOCA time.  sign of DOCA is irrelevant here.
                double hitt0 = tdrift - doca/_vdrift;
		//		std::cout<<"t0fit "<<hitt0-myfit._t00.t0()<<" "<<doca<<std::endl;
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
	    //	    cout<<"t0fitselect "<<niter<<" "<<t0<<" "<<t0sig<<" "<<nactive<<endl;
	  }
	}
// reset t0
//	t0err=0;
        myfit._t0.setT0(t0,t0err);
	if(myfit._krep != 0)myfit._krep->setT0(myfit._t0);
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
  KalFitI::weedHits(TrkKalFit& myfit) {
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
      worstHot->setUsability(1);
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

  bool
  KalFitI::unweedHits(TrkKalFit& myfit) {
    // Loop over HoTs and find HoT with largest contribution to chi2.  If this value
    // is greater than some cut value, deactivate that HoT and reFit
    bool retval(false);
    double worst = 1.e12;
    TrkHitOnTrk* worstHot = 0;
    TrkHotList* hots = myfit._krep->hotList();
    for (TrkHotList::nc_hot_iterator iHot = hots->begin(); iHot != hots->end(); ++iHot) {
      if (!iHot->isActive()) {
        double resid, residErr;
        if(iHot->resid(resid, residErr, true)){
          double value = fabs(resid/residErr);
	  if(_debug>1) {
	    TrkHitOnTrk *hit=iHot.get();
	    TrkStrawHit *shit=(TrkStrawHit*)hit;
	    cout<<shit->index()<<" unwind chi2 "<<value<<" max "<<_maxhitchi<<endl;
	  }
          if (value < _maxhitchi && value < worst) {
            worst = value;
            worstHot = iHot.get();
          }
        }
      }
    }
    if(0 != worstHot){
      retval = true;
      worstHot->setActivity(true);
      worstHot->setUsability(3);
      myfit.fit();
      myfit._krep->addHistory(myfit._fit, "HitUnWeed");
      // Recursively iterate
      myfit._nweediter++;
      if (myfit._fit.success() && myfit._nweediter < _maxweed ) {
        retval |= unweedHits(myfit);
      }
    }
    return retval;
  }

  const BField*
  KalFitI::bField() {
    if(_bfield == 0){
// create a wrapper around the mu2e nominal DS field
      GeomHandle<BFieldConfig> bfconf;
      _bfield=new BFieldFixed(bfconf->getDSUniformValue());
      assert(_bfield != 0);
    }
    return _bfield;
  }

  double KalFitI::findZFltlen(const TrkKalFit& myfit,double zval) {
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

  void KalFitI::reActivateHitsbyChi2(TrkDef& /*mytrk*/,TrkKalFit& myfit){
    unsigned niter=0;
    bool changed=true;
    while(changed && myfit._fit.success() && niter < _kalcon->maxIterations()){
      changed = false;
      // re add by chi2
      if(_weedhits){
	myfit._nweediter = 0;
	changed |= unweedHits(myfit);
	std::cout<<"afterunweedchi2 "<<myfit._krep->chisquared(trkIn)<<" t0 "<<myfit._t0.t0()<<" ndof "<<myfit._krep->nDof()<<std::endl;
      }
      niter++;
    }
  }

  void KalFitI::reActivateHitsbyTurn(TrkDef& mytrk,TrkKalFit& myfit){
    double midflt = 0.5*(myfit._krep->lowFitRange() + myfit._krep->hiFitRange());
    double locflt;
    const HelixTraj* shelix = dynamic_cast<const HelixTraj*>(myfit._krep->localTrajectory(midflt,locflt));
    mytrk.setHelix(*shelix);
    mytrk.setTraj(&myfit._krep->pieceTraj());
    int nfixed=0;
    for(std::vector<TrkStrawHit*>::iterator ihit=myfit._hits.begin();ihit!=myfit._hits.end();ihit++){
      TrkStrawHit* trkhit = *ihit;
      if (!trkhit->isActive()) {
	if(fixHitTurn(mytrk,trkhit))
	  nfixed++;
      }
    }
    std::cout<<"number of fixed hits by turn "<<nfixed<<std::endl;
    if(nfixed){
      myfit._krep->resetFit();
      fitTrack(myfit);
    }
  }

}
