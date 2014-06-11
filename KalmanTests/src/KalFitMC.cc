//
// MC functions associated with KalFit
// $Id: KalFitMC.cc,v 1.59 2014/06/11 00:20:14 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/06/11 00:20:14 $
//
//geometry
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "BFieldGeom/inc/BFieldConfig.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
// services
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "art/Framework/Services/Optional/TFileService.h"
// data
#include "art/Framework/Principal/Event.h"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/VirtualDetectorId.hh"
// Utilities
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
// tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
// BaBar
#include "BaBar/BaBar.hh"
#include "KalmanTests/inc/TrkDef.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTests/inc/KalFitMC.hh"
#include "TrkBase/TrkHelixUtils.hh"
#include "TrkBase/TrkHotList.hh"
#include "TrkBase/TrkPoca.hh"
#include "TrkBase/TrkMomCalculator.hh"
#include "BField/BFieldFixed.hh"
#include "KalmanTrack/KalHit.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// root 
#include "TMath.h"
#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TGMsgBox.h"
#include "TTree.h"
// C++
#include <iostream>
#include <string>
#include <functional>

using namespace std; 
 
namespace mu2e 
{
  // comparison functor for ordering step points
  struct timecomp : public binary_function<MCStepItr,MCStepItr, bool> {
    bool operator()(MCStepItr x,MCStepItr y) { return x->time() < y->time(); }
  };
  struct devicecomp : public binary_function<const TrkStrawHit*, const TrkStrawHit*, bool> {
    bool operator()(const TrkStrawHit* x, const TrkStrawHit* y) { 
// predicate on device first, as fltlen might be ambiguous for inactive hots
      if(x->straw().id().getDevice() == y->straw().id().getDevice())
	return x->fltLen() < y->fltLen();
      else
	return(x->straw().id().getDevice() < y->straw().id().getDevice());
    }
  };

  struct spcount {
    spcount() : _count(0) {}
    spcount(art::Ptr<SimParticle> const& spp) : _spp(spp), _count(1) {}
    void append(art::Ptr<SimParticle> const& sp) { if(sp == _spp)++_count; }
    bool operator ==(art::Ptr<SimParticle> const& sp) const { return _spp == sp; }
    art::Ptr<SimParticle> _spp;
    unsigned _count;
  };

  MCHitSum::MCHitSum(StepPointMC const& mchit,art::Ptr<SimParticle>& spp) :
    _esum(mchit.eDep()),_count(1),
  _spp(spp),_pdgid(0),_gid(-1),_pid(-1),_sid(mchit.strawIndex()),_mcsid(mchit.strawIndex()),
  _t0(mchit.time()),_time(mchit.time()),
  _pos(mchit.position()),
  _mom(mchit.momentum()){
    if(spp.isNonnull() ){
      _pdgid = spp->pdgId();
      _pid = spp->creationCode();
      if( spp->genParticle().isNonnull())
	_gid = spp->genParticle()->generatorId().id();
    }
  }
  MCHitSum::MCHitSum(StrawDigiMC const& mcdigi) : _esum(mcdigi.energySum()),
  _count(mcdigi.stepPointMCs().size()),
  _pdgid(0),_gid(-1),_pid(-1),_sid(mcdigi.strawIndex()),_mcsid(0),_time(-1000.0)
  {
  // primary end is end '0'
    StrawDigi::TDCChannel itdc = StrawDigi::zero;
    if(mcdigi.hasTDC(itdc)){
      _mcsid = mcdigi.stepPointMCs()[itdc]->strawIndex();
      _pos = mcdigi.stepPointMC(itdc)->position();
      _mom = mcdigi.stepPointMC(itdc)->momentum();
      _t0 = mcdigi.clusterPosition(itdc).t();
      _time = mcdigi.wireEndTime(itdc);
      _pos = mcdigi.clusterPosition(itdc);
      _spp = mcdigi.stepPointMC(itdc)->simParticle();
      if(_spp.isNonnull() ){
	_pdgid = _spp->pdgId();
	_pid = _spp->creationCode();
	if( _spp->genParticle().isNonnull())
	  _gid = _spp->genParticle()->generatorId().id();
      }
    }
  }
  void MCHitSum::append(StepPointMC const& mchit)  {
    double eold = _esum;
    _esum += mchit.eDep();
    _time = _time*eold/_esum + mchit.time()*mchit.eDep()/_esum;
    _pos = _pos*(eold/_esum) + mchit.position()*(mchit.eDep()/_esum);
    _mom = _mom*(eold/_esum) + mchit.momentum()*(mchit.eDep()/_esum);
    _count++;
  }

  KalFitMC::~KalFitMC(){}
  
  KalFitMC::KalFitMC(fhicl::ParameterSet const& pset) :
    _mcptrlabel(pset.get<std::string>("MCPtrLabel","makeSH")),
    _mcstepslabel(pset.get<std::string>("MCStepsLabel","g4run")),
    _simpartslabel(pset.get<std::string>("SimParticleLabel","g4run")),
    _mcdigislabel(pset.get<std::string>("StrawHitMCLabel","makeSH")),
    _strawhitslabel(pset.get<std::string>("strawHitsLabel","makeSH")),
    _mintrkmom(pset.get<double>("minTrkMom",60.0)),
    _mct0err(pset.get<double>("mcT0Err",0.1)),
    _mcambig(pset.get<bool>("mcAmbiguity",false)),
    _debug(pset.get<int>("debugLevel",0)),
    _diag(pset.get<int>("diagLevel",1)),
    _minnhits(pset.get<unsigned>("minNHits",10)),
    _maxnhits(pset.get<unsigned>("maxNHits",120)),
    _maxarcgap(pset.get<int>("MaxArcGap",2)),
    _purehits(pset.get<bool>("pureHits",false)),
    _trkdiag(0),_hitdiag(0)
  {
// define the ids of the virtual detectors
    _midvids.push_back(VirtualDetectorId::TT_Mid);
    _midvids.push_back(VirtualDetectorId::TT_MidInner);
    _entvids.push_back(VirtualDetectorId::TT_FrontHollow);
    _entvids.push_back(VirtualDetectorId::TT_FrontPA); 
    _xitvids.push_back(VirtualDetectorId::TT_Back);
  }

// Find MC truth for the given particle entering a given detector(s).  Unfortunately the same logical detector can map
// onto several actual detectors, as they are not allowed to cross volume boundaries, so we have to check against
// several ids.  The track may also pass through this detector more than once, so we return a vector, sorted by time.
  void 
  KalFitMC::findMCSteps(StepPointMCCollection const* mcsteps, cet::map_vector_key const& trkid, std::vector<int> const& vids,
  std::vector<MCStepItr>& steps) {
    steps.clear();
    if(mcsteps != 0){
      // Loop over the step points, and find the one corresponding to the given detector
      for( MCStepItr imcs =mcsteps->begin();imcs!= mcsteps->end();imcs++){
	if(vids.size() == 0 ||  (imcs->trackId() == trkid && std::find(vids.begin(),vids.end(),imcs->volumeId()) != vids.end())){
	  steps.push_back(imcs);
	}
      }
      // sort these in time
      std::sort(steps.begin(),steps.end(),timecomp());
    }
  }

// define seed helix, t0, and hits coming from a given particle using MC truth.  Note that the input
// trkdef object must reference a valid straw hit collection
  bool
  KalFitMC::trkFromMC(cet::map_vector_key const& trkid,TrkDef& mytrk) {
    if(_mcdata._mcsteps == 0)return false;
// preset to failure
    bool retval(false);
    GlobalConstantsHandle<ParticleDataTable> pdt;
    const Tracker& tracker = getTrackerOrThrow();
    unsigned nstraws = mytrk.strawHitCollection()->size();
    if(nstraws >= _minnhits ){
      GeomHandle<DetectorSystem> det;
      GeomHandle<VirtualDetector> vdg;
      std::vector<hitIndex> indices;
// find the mcstep at the middle of the detector
      std::vector<MCStepItr> steps;
      findMCSteps(_mcdata._mcvdsteps,trkid,_midvids,steps);
      if(steps.size() > 0 && vdg->exist(steps[0]->volumeId())&& steps[0]->momentum().mag() > _mintrkmom){
// take the first point
        MCStepItr imcs = steps[0];
        double t0 = imcs->time();
        double charge = pdt->particle(imcs->simParticle()->pdgId()).ref().charge();
        CLHEP::Hep3Vector mom = imcs->momentum();
// need to transform into the tracker coordinate system
        CLHEP::Hep3Vector pos = det->toDetector(imcs->position());
        if(_debug > 1)std::cout << "Defining track at virtual detector id= " << imcs->volumeId()
          << " name " << VirtualDetectorId(imcs->volumeId()).name()
          << " position = " << pos
          << " momentum = " << mom
          << " time = " << t0 << std::endl;
// find the indices of the true primary particle straw hits
        for(size_t istraw=0;istraw<nstraws;istraw++){
          PtrStepPointMCVector const& mcptr(_mcdata._mchitptr->at(istraw));
	  StrawHit const& sh = mytrk.strawHitCollection()->at(istraw);
	  Straw const& straw = tracker.getStraw(sh.strawIndex());
	  int ambig(0);
	  if(_mcambig){
// compute the hit ambiguity using MC truth, using the earliest energy deposit
	    Hep3Vector dir = mcptr[0]->momentum().unit();
	    Hep3Vector start = mcptr[0]->position();
	    Hep3Vector end = start + dir*mcptr[0]->stepLength();
	    Hep3Vector mid = 0.5*(start+end);
	    Hep3Vector mcsep = straw.getMidPoint() - mid;
	    Hep3Vector mcperp = (dir.cross(straw.getDirection())).unit();
	    double dperp = mcperp.dot(mcsep);
	    ambig = dperp > 0 ? 1 : -1; // follow TrkPoca convention
	  }
	  unsigned nprimary(0);
          for( size_t j=0; j<mcptr.size(); ++j ) {
            StepPointMC const& mchit = *mcptr[j];
 // not sure if this works with bkgs merged FIXME!!!!  we also need to treat energy from daughters different
 // from energy from completely different particles
            if( mchit.trackId() == trkid )nprimary++;
          }
// decide if we want all hits with any contribution from this particle, or only pure hits from this particle.
          if( nprimary > 0 && (nprimary == mcptr.size() || (!_purehits) ) )indices.push_back(hitIndex(istraw,ambig));
        }
        if(indices.size() >= _minnhits && indices.size() <= _maxnhits){
// nominal magnetic field.
          HepVector parvec(5,0);
// Babar interface still uses HepPoint: FIXME!!!!
// use the z component of th enominal field to define the initial helix parameters.  Off-axis terms are
// ignored anyways by this interface
          double hflt(0.0);
	  GeomHandle<BFieldManager> bfmgr;
	  GeomHandle<DetectorSystem> det;
	  CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(Hep3Vector(0.0,0.0,0.0));
	  double bz = bfmgr->getBField(vpoint_mu2e).z();
	  TrkHelixUtils::helixFromMom( parvec, hflt, 
            HepPoint(pos.x(),pos.y(),pos.z()),
            mom,charge,bz);
  // dummy covariance matrix; this should be set according to measured values, FIXME!!!!!
          HepSymMatrix dummy(5,1); 
          dummy(1,1)=1.; dummy(2,2)=0.1*0.1;dummy(3,3)=1e-2*1e-2;
          dummy(4,4)=1.; dummy(5,5)=0.1*0.1;
	  TrkFitDirection::FitDirection fdir = parvec[4]>0 ?
	  TrkFitDirection::downstream : TrkFitDirection::upstream;
	  TrkParticle::type ptype = (TrkParticle::type)imcs->simParticle()->pdgId();
          mytrk = TrkDef(mytrk.strawHitCollection(),indices,parvec,dummy,TrkParticle(ptype),TrkFitDirection(fdir));
	  mytrk.setT0(TrkT0(t0,_mct0err));
          retval = true;
        }
      }
    }
    return retval;
  }
 
  void KalFitMC::hitDiag(const TrkStrawHit* strawhit) {
    if(_hitdiag == 0)createHitDiag();
// straw info
    _edep = strawhit->strawHit().energyDep();
//    CLHEP::Hep3Vector dmid = strawhit->wirePosition()-strawhit->straw().getMidPoint(); 
// trkhit  info
    _rdrift = strawhit->driftRadius();
    _rdrifterr = strawhit->driftRadiusErr();
    _amb = strawhit->ambig();
    CLHEP::Hep3Vector shpos;
    strawhit->hitPosition(shpos);
    _shpos = shpos;
//    _shposs.push_back(shpos);
    _dmid = strawhit->timeDiffDist();
    _dmiderr = strawhit->timeDiffDistErr();
// residual info
    double resid,residerr;
    if(strawhit->resid(resid,residerr,true)){
      _resid = resid;
      _residerr = residerr;
    } else {
      _resid = _residerr = -100.0;
    }
    _hitt0 = strawhit->hitT0()._t0;
    _hitt0err = strawhit->hitT0()._t0err;
    _hflt = strawhit->hitLen();
    _trkflt = strawhit->fltLen();
    _active = strawhit->isActive();
    _use = strawhit->usability();
// mc information
    unsigned istraw = strawhit->index();

    CLHEP::Hep3Vector mcpos;
    PtrStepPointMCVector const& mcptr(_mcdata._mchitptr->at(istraw));
    _nmcsteps = mcptr.size();
    if(_nmcsteps > 0){
      const Tracker& tracker = getTrackerOrThrow();
      Straw const& straw = tracker.getStraw(mcptr[0]->strawIndex());
      TwoLinePCA pca( straw.getMidPoint(), straw.getDirection(), 
	  mcptr[0]->position(), mcptr[0]->momentum().unit() );
      _mcrdrift = pca.dca();
      _mcdmid = pca.s1();
    }
// sums
    double esum(0.0);
    double mct0(0.0);
    _pdist = 0.0;
    _pperp = 0.0;
    _pmom = 0.0;
    for( size_t imc=0; imc< _nmcsteps; ++imc ) {
      StepPointMC const& mchit = *mcptr[imc];
      if( mchit.trackId().asInt() == 1 ){
        double edep = mchit.eDep();
        esum += edep;
        mcpos += mchit.position()*edep;
        mct0 += mchit.time()*edep;
        CLHEP::Hep3Vector dprod = mchit.position()-mchit.simParticle()->startPosition();
        _pdist += dprod.mag()*edep;
        _pperp += dprod.perp(Hep3Vector(0.0,0.0,1.0))*edep;
        _pmom = mchit.momentum().mag()*edep;
      }
    }
    if(esum > 0.0){
      mcpos /= esum;
      mct0 /= esum;
      _pdist /= esum;
      _pperp /= esum;
      _pmom /= esum;
    }
    _mchpos = mcpos;
    _mchitt0 = mct0;
    _hitdiag->Fill();
  }
  
  void
  KalFitMC::fillHitsVector(const KalRep* krep,std::vector<const TrkStrawHit*> & hits) {
    hits.clear();
    if(krep != 0) {
  // extract the hits from the kalrep and perform diagnstics
      const TrkHotList* hots = krep->hotList();
      hits.reserve(hots->nHit());
      for(TrkHotList::hot_iterator ihot=hots->begin();ihot != hots->end();++ihot){
	const TrkStrawHit* hit = dynamic_cast<const TrkStrawHit*>(ihot.get());
	if(hit != 0)hits.push_back(hit);
      }
      std::sort(hits.begin(),hits.end(),devicecomp());
    }
  }

  void
  KalFitMC::kalDiag(const KalRep* krep,bool fill) {
    GeomHandle<VirtualDetector> vdg;
    GeomHandle<DetectorSystem> det;

    if(krep != 0) {
      if(_trkdiag == 0)createTrkDiag();
// find the associated MC partcle
      art::Ptr<SimParticle> spp;
      findMCTrk(krep,spp);
// mc track patermeter info for the particle which generated most of the hits
      if(spp.isNonnull() ){
	SimParticle const& sp = *spp;
	mcTrkInfo(sp);
	if( sp.genParticle().isNonnull())
	  _mcgenid = sp.genParticle()->generatorId().id();
	else
	  _mcgenid = 0;
	_mcpdgid = sp.pdgId();
	_mcproc = sp.creationCode();
      } else {
	_mcpdgid = 0;
	_mcgenid = -1;
	_mcproc = -1;
      }
// no information on iterations either!
      _nt0iter = _nweediter = -1;
      if(_diag > 1){
	std::vector<const TrkStrawHit*> hits;
	fillHitsVector(krep,hits);	
	hitsDiag(hits);
      }
      if(krep->fitCurrent()){
	_t0 = krep->t0().t0();
	_t0err = krep->t0().t0Err();
	_fitstatus = krep->fitStatus().success();
	_nhits = krep->hotList()->nHit();
	_niter = krep->iterations();
	_ndof = krep->nDof();
	_nactive = krep->nActive();
	_chisq = krep->chisq();
	_fitcon = krep->chisqConsistency().significanceLevel();
	_radlen = krep->radiationFraction();
	_nsites = krep->siteList().size();
	_firstflt = krep->firstHit()->globalLength();
	_lastflt = krep->lastHit()->globalLength();
	// get the fit at the first hit
	CLHEP::Hep3Vector entpos = det->toDetector(vdg->getGlobal(VirtualDetectorId::TT_FrontPA));
	double zent = entpos.z();
	double firsthitfltlen = krep->firstHit()->kalHit()->hitOnTrack()->fltLen() - 10;
	double lasthitfltlen = krep->lastHit()->kalHit()->hitOnTrack()->fltLen() - 10;
	double entlen = std::min(firsthitfltlen,lasthitfltlen);
	TrkHelixUtils::findZFltlen(krep->traj(),zent,entlen,0.1); 
	double loclen(0.0);
	const TrkSimpTraj* ltraj = krep->localTrajectory(entlen,loclen);
	_fitpar = helixpar(ltraj->parameters()->parameter());
	_fiterr = helixpar(ltraj->parameters()->covariance());
	CLHEP::Hep3Vector fitmom = krep->momentum(entlen);
	BbrVectorErr momerr = krep->momentumErr(entlen);
	_fitmom = fitmom.mag();
	Hep3Vector momdir = fitmom.unit();
	HepVector momvec(3);
	for(int icor=0;icor<3;icor++)
	  momvec[icor] = momdir[icor];
	_fitmomerr = sqrt(momerr.covMatrix().similarity(momvec));
	CLHEP::Hep3Vector seedmom = TrkMomCalculator::vecMom(*(krep->seed()),krep->kalContext().bField(),0.0);
	_seedmom = seedmom.mag();
      } else {
	_fitstatus = -krep->fitStatus().failure();
	_nhits = -1;
	_fitmom = -1.0;
	_seedmom = -1.0;
      }
    } else {
      _fitstatus = -1000;
      _nhits = -1;
      _fitmom = -1.0;
      _seedmom = -1.0;
      _mcpdgid = 0;
      _mcgenid = -1;
      _mcproc = -1;
// use the primary particle (if it exists)
      static cet::map_vector_key trkid(1);
      if(_mcdata._simparts != 0){
	for ( SimParticleCollection::const_iterator isp = _mcdata._simparts->begin();
	    isp != _mcdata._simparts->end(); ++isp ){
	  if(isp->second.id() == trkid){
	    SimParticle const& sp = isp->second;
	    _mcgenid = sp.genParticle()->generatorId().id();
	    _mcpdgid = sp.pdgId();
	    _mcproc = sp.creationCode();
	    fillMCTrkInfo(sp,_mcinfo);
	    break;
	  }
	}
      }
    }
    if(fill)_trkdiag->Fill(); 
  }


  void
  KalFitMC::findMCTrk(const KalRep* krep,art::Ptr<SimParticle>& spp) {
// find the SimParticle which contributed most of the hits
   std::vector<spcount> sct;
// get the straw hits from the track
    std::vector<const TrkStrawHit*> hits;
    fillHitsVector(krep,hits);	
// loop over the hits and find the associated steppoints
    for(size_t itsh=0;itsh<hits.size();++itsh){
      const TrkStrawHit* tsh = hits[itsh];
      PtrStepPointMCVector const& mcptr(_mcdata._mchitptr->at(tsh->index()));
      std::vector<MCHitSum> mcinfo;
      fillMCHitSum(mcptr,mcinfo);
      if(mcinfo.size() > 0){
	bool found(false);
	for(size_t isp=0;isp<sct.size();++isp){
	  if(sct[isp] == mcinfo[0]._spp){
	    sct[isp].append(mcinfo[0]._spp);
	    found = true;
	    break;
	  }
	}
	if(!found)sct.push_back(spcount(mcinfo[0]._spp));
      }
    }
    unsigned imax(0);
    for(size_t isp=0;isp<sct.size();++isp){
      if(sct[isp]._count>imax){
	imax = sct[isp]._count;
	spp = sct[isp]._spp;
      }
    }
  }
  
  void KalFitMC::hitsDiag(std::vector<const TrkStrawHit*> const& hits) {
    _tshinfo.clear();
    _ncactive = 0;
 // loop over hits
    for(size_t itsh=0;itsh<hits.size();++itsh){
      const TrkStrawHit* tsh = hits[itsh];
      if(tsh != 0){
        TrkStrawHitInfo tshinfo;
        tshinfo._active = tsh->isActive();
        tshinfo._usable = tsh->usability();
	tshinfo._device = tsh->straw().id().getDevice();
	tshinfo._sector = tsh->straw().id().getSector();
	tshinfo._layer = tsh->straw().id().getLayer();
	tshinfo._straw = tsh->straw().id().getStraw();
	tshinfo._edep = tsh->strawHit().energyDep();
	static HepPoint origin(0.0,0.0,0.0);
	CLHEP::Hep3Vector hpos = tsh->hitTraj()->position(tsh->hitLen()) - origin;
	tshinfo._z = hpos.z();
	tshinfo._phi = hpos.phi();
	tshinfo._rho = hpos.perp();
	double resid,residerr;
	if(tsh->resid(resid,residerr,true)){
	  tshinfo._resid = resid;
	  tshinfo._residerr = residerr;
	} else {
	  tshinfo._resid = tshinfo._residerr = -1000.;
	}
	tshinfo._rdrift = tsh->driftRadius();
	tshinfo._rdrifterr = tsh->driftRadiusErr();
	double rstraw = tsh->straw().getRadius();
	tshinfo._dx = sqrt(std::max(0.0,rstraw*rstraw-tshinfo._rdrift*tshinfo._rdrift));
	tshinfo._trklen = tsh->fltLen();
	tshinfo._hlen = tsh->hitLen();
	tshinfo._t0 = tsh->hitT0()._t0;
	// include signal propagation time correction
	tshinfo._ht = tsh->time()-tsh->signalTime();
	tshinfo._tddist = tsh->timeDiffDist();
	tshinfo._tdderr = tsh->timeDiffDistErr();
	tshinfo._ambig = tsh->ambig();
	if(tsh->poca() != 0)
	  tshinfo._doca = tsh->poca()->doca();
	else
	  tshinfo._doca = -100.0;
	tshinfo._exerr = tsh->extErr();
	tshinfo._penerr = tsh->penaltyErr();
	tshinfo._t0err = tsh->t0Err()/tsh->driftVelocity();
	// MC information	
	PtrStepPointMCVector const& mcptr(_mcdata._mchitptr->at(tsh->index()));
	tshinfo._mcn = mcptr.size();
	std::vector<MCHitSum> mcsum;
	if(_mcdata._mcdigis != 0)
	  mcsum.push_back(_mcdata._mcdigis->at(tsh->index()));
	else if(_mcdata._mcsteps != 0)
	  fillMCHitSum(mcptr,mcsum);
	tshinfo._mcnunique = mcsum.size();
	if(mcsum.size() > 0){
	  tshinfo._mcppdg = mcsum[0]._pdgid;
	  tshinfo._mcpgen = mcsum[0]._gid;
	  tshinfo._mcpproc = mcsum[0]._pid;
// first hit is the one that set t0
	  tshinfo._mct0 = mcsum[0]._t0;
	  tshinfo._mcht = mcsum[0]._time;
	  art::Ptr<SimParticle> sp = mcsum[0]._spp;
	  tshinfo._mcpdg = sp->pdgId();
	  tshinfo._mcproc = sp->creationCode();
	  tshinfo._mcedep = mcsum[0]._esum;
	  tshinfo._mcgen = -1;
	  if(sp->genParticle().isNonnull())
	    tshinfo._mcgen = sp->genParticle()->generatorId().id();
// find the step midpoint
	  Hep3Vector mcsep = mcsum[0]._pos-tsh->straw().getMidPoint();
	  Hep3Vector dir = mcsum[0]._mom.unit();
	  Hep3Vector mcperp = (dir.cross(tsh->straw().getDirection())).unit();
	  double dperp = mcperp.dot(mcsep);
	  tshinfo._mcdist = fabs(dperp);
	  tshinfo._mcambig = dperp > 0 ? -1 : 1; // follow TrkPoca convention
	  tshinfo._mclen = mcsep.dot(tsh->straw().getDirection());
	  tshinfo._xtalk = mcsum[0]._sid != tsh->strawHit().strawIndex();
	}
	_tshinfo.push_back(tshinfo);
// count active conversion hits
	if(tshinfo._mcgen==2&&tsh->isActive())++_ncactive;
      }
    }
  }

  void KalFitMC::mcTrkInfo(SimParticle const&  sp) {
    GeomHandle<VirtualDetector> vdg;
    GeomHandle<DetectorSystem> det;
// find the mc info at the entrance to the detector
    fillMCTrkInfo(sp,_mcinfo);
    // find MC info at tracker
    cet::map_vector_key trkid = sp.id();
    std::vector<MCStepItr> entsteps,midsteps,xitsteps;
    findMCSteps(_mcdata._mcvdsteps,trkid,_entvids,entsteps);
    if(entsteps.size() > 0 && vdg->exist(entsteps[0]->volumeId()))
      fillMCTrkInfo(entsteps.front(),_mcentinfo);
    findMCSteps(_mcdata._mcvdsteps,trkid,_midvids,midsteps);
    if(midsteps.size() > 0 && vdg->exist(midsteps[0]->volumeId()))
      fillMCTrkInfo(midsteps.front(),_mcmidinfo);
    findMCSteps(_mcdata._mcvdsteps,trkid,_xitvids,xitsteps);
    if(xitsteps.size() > 0 && vdg->exist(xitsteps[0]->volumeId()))
      fillMCTrkInfo(xitsteps.front(),_mcxitinfo);
// find MC info about brems in the tracker
// should get these numbers from a service, FIXME!!!!
    static double _trkminz(-1501.0);
    static double _trkmaxz(1501.0);
    _bremsesum = 0.0;
    _bremsemax = 0.0;
    _bremsz = _trkminz;
    if(_mcdata._simparts != 0){
      for ( SimParticleCollection::const_iterator isp = _mcdata._simparts->begin();
	  isp != _mcdata._simparts->end(); ++isp ){
	SimParticle const& sp = isp->second;
	// find photons with parent = the conversion electron created by brems
	if(sp.hasParent() && sp.parent()->id() == trkid && sp.pdgId() == PDGCode::gamma  && sp.creationCode() == ProcessCode::eBrem){
	  CLHEP::Hep3Vector pos = det->toDetector(sp.startPosition());
	  if(pos.z() > _trkminz && pos.z() < _trkmaxz){
	    double de = sp.startMomentum().e();
	    if(de > _bremsemax){
	      _bremsemax = de;
	      _bremsz = pos.z();
	    }
	    _bremsesum += de;
	  }
	}
      }
    }
  }

  void KalFitMC::fillMCHitSummary(){
    size_t nsh = _mcdata._mchitptr->size();
    _mchitsums.reserve(nsh);
    for(size_t ish=0;ish<nsh;++ish){
//      PtrStepPointMCVector const& mcptr(_mcdata._mchitptr->at(ish));
      std::vector<MCHitSum> mcsum;
//      fillMCHitSum(mcptr,mcsum);
      if(mcData()._mcdigis != 0)
	mcsum.push_back(MCHitSum(mcData()._mcdigis->at(ish)));
      _mchitsums.push_back(mcsum);
    }
  }

// Summarize an associated set of StepPointMCs from a StrawHit according to their parents.  This assignes
// daughter contributions to their parents
  void KalFitMC::fillMCHitSum(PtrStepPointMCVector const& mcptr,std::vector<MCHitSum>& summary ){
// first, create a map from daughters to mothers
    std::map<SPPtr,SPPtr> mdmap;
    findRelatives(mcptr,mdmap);
// Loop over the step points, and fill the summary vector
    for( size_t imc=0; imc< mcptr.size(); ++imc ) {
// find the primary parent, and create a pair for this step point's energy deposition
      art::Ptr<SimParticle> sp = mcptr[imc]->simParticle();
// find it's parent
      art::Ptr<SimParticle> spp = mdmap[sp];
// create the summary
      MCHitSum tsum(*mcptr[imc],spp);
// Add this energy to this particle, or create the entry if this is the first time this particle is seen
      std::vector<MCHitSum>::iterator ifnd = std::find(summary.begin(),summary.end(),tsum);
      if(ifnd == summary.end())
        summary.push_back(tsum);
      else
        ifnd->append(*mcptr[imc]);
    }
// sort this according to deposited energy
    std::sort(summary.begin(),summary.end(),MCHitSum::ecomp());
  }

// map daughters onto parents within the context of an associated set of StepPointMCs (like from a StrawHit).
  void KalFitMC::findRelatives(PtrStepPointMCVector const& mcptr,std::map<SPPtr,SPPtr>& mdmap){
    // loop over steps
    for( size_t imc=0; imc< mcptr.size(); ++imc ) {
      art::Ptr<SimParticle> spp = mcptr[imc]->simParticle();		  
      if(mdmap.find(spp) == mdmap.end()){
    // start with self-reference
        mdmap[spp] = spp;
      }
    }
    // loop over the steps again, trying to map particles to their parents
    for( size_t imc=0; imc< mcptr.size(); ++imc ) {
      art::Ptr<SimParticle> spp = mcptr[imc]->simParticle();		  
      if(mdmap[spp] == spp && !spp.isNull()){
    // move up through the genealogy till we find the highest-rank parent directly contributing
    // to this strawhit
        art::Ptr<SimParticle> sppp = spp->parent();
        while(!sppp.isNull()){
          std::map<SPPtr,SPPtr>::iterator ifnd = mdmap.find(sppp);
          if(ifnd != mdmap.end())
    // repoint this particle to its highest-level contributing parent
            mdmap[spp] = sppp;
    // move on to the parent's parent, as there can be gaps in the direct contribution
          sppp = sppp->parent();
        }
      }
    }
  }
 
  TTree*
  KalFitMC::createTrkDiag(){
    art::ServiceHandle<art::TFileService> tfs;
    _trkdiag=tfs->make<TTree>("trkdiag","trk diagnostics");
    _trkdiag->Branch("fitstatus",&_fitstatus,"fitstatus/I");
    _trkdiag->Branch("t0",&_t0,"t0/F");
    _trkdiag->Branch("t0err",&_t0err,"t0err/F");
    _trkdiag->Branch("nhits",&_nhits,"nhits/I");
    _trkdiag->Branch("ndof",&_ndof,"ndof/I");
    _trkdiag->Branch("niter",&_niter,"niter/I");
    _trkdiag->Branch("nt0iter",&_nt0iter,"nt0iter/I");
    _trkdiag->Branch("nweediter",&_nweediter,"nweediter/I");
    _trkdiag->Branch("nactive",&_nactive,"nactive/I");
    _trkdiag->Branch("ncactive",&_ncactive,"ncactive/I");
    _trkdiag->Branch("narcs",&_narcs,"narcs/I");
    _trkdiag->Branch("nchits",&_nchits,"nchits/I");
    _trkdiag->Branch("ncgood",&_ncgood,"ncgood/I");
    _trkdiag->Branch("chisq",&_chisq,"chisq/F");
    _trkdiag->Branch("fitcon",&_fitcon,"fitcon/F");
    _trkdiag->Branch("radlen",&_radlen,"radlen/F");
    _trkdiag->Branch("firstflt",&_firstflt,"firstflt/F");
    _trkdiag->Branch("lastflt",&_lastflt,"lastflt/F");
    _trkdiag->Branch("nsites",&_nsites,"nsites/I");
    _trkdiag->Branch("fitmom",&_fitmom,"fitmom/F");
    _trkdiag->Branch("fitmomerr",&_fitmomerr,"fitmomerr/F");
    _trkdiag->Branch("seedmom",&_seedmom,"seedmom/F");
    _trkdiag->Branch("fitpar",&_fitpar,"d0/F:p0/F:om/F:z0/F:td/F");
    _trkdiag->Branch("fiterr",&_fiterr,"d0err/F:p0err/F:omerr/F:z0err/F:tderr/F");
// basic MC info
    _trkdiag->Branch("mcpdgid",&_mcpdgid,"mcpdgid/I");
    _trkdiag->Branch("mcgenid",&_mcgenid,"mcgenid/I");
    _trkdiag->Branch("mcproc",&_mcproc,"mcproc/I");
// mc info at production and several spots in the tracker
    _trkdiag->Branch("mcinfo",&_mcinfo,"mcpdgid/I:mct0/F:mcmom/F:mcx/F:mcy/F:mcz/F:mcd0/F:mcp0/F:mcom/F:mcz0/F:mctd/F");
    _trkdiag->Branch("mcentinfo",&_mcentinfo,"mcentpdgid/I:mcentt0/F:mcentmom/F:mcentx/F:mcenty/F:mcentz/F:mcentd0/F:mcentp0/F:mcentom/F:mcentz0/F:mcenttd/F");
    _trkdiag->Branch("mcmidinfo",&_mcmidinfo,"mcmidpdgid/I:mcmidt0/F:mcmidmom/F:mcmidx/F:mcmidy/F:mcmidz/F:mcmidd0/F:mcmidp0/F:mcmidom/F:mcmidz0/F:mcmidtd/F");
    _trkdiag->Branch("mcxitinfo",&_mcxitinfo,"mcxitpdgid/I:mcxitt0/F:mcxitmom/F:mcxitx/F:mcxity/F:mcxitz/F:mcxitd0/F:mcxitp0/F:mcxitom/F:mcxitz0/F:mcxittd/F");
// info about energy loss in the tracker
    _trkdiag->Branch("bremsesum",&_bremsesum,"bremsesum/F");
    _trkdiag->Branch("bremsemax",&_bremsemax,"bremsemax/F");
    _trkdiag->Branch("bremsz",&_bremsz,"bremsz/F");
// track hit and arc info    
    if(_diag > 1)
      _trkdiag->Branch("tshinfo",&_tshinfo);
    return _trkdiag;
  }

  TTree*
  KalFitMC::createHitDiag(){
    art::ServiceHandle<art::TFileService> tfs;
    _hitdiag=tfs->make<TTree>("hitdiag","hit diagnostics");      
    _hitdiag->Branch("shpos",&_shpos,"x/F:y/F:z/F");
    _hitdiag->Branch("rdrift",&_rdrift,"rdrift/F");
    _hitdiag->Branch("rdrifterr",&_rdrifterr,"rdrifterr/F");
    _hitdiag->Branch("amb",&_amb,"amb/I");
    _hitdiag->Branch("hitt0",&_hitt0,"hitt0/F");
    _hitdiag->Branch("hitt0err",&_hitt0err,"hitt0err/F");
    _hitdiag->Branch("dmid",&_dmid,"dmid/F");
    _hitdiag->Branch("dmiderr",&_dmiderr,"dmiderr/F");
    _hitdiag->Branch("resid",&_resid,"resid/F");
    _hitdiag->Branch("residerr",&_residerr,"residerr/F");
    _hitdiag->Branch("hflt",&_hflt,"hflt/F");
    _hitdiag->Branch("trkflt",&_trkflt,"trkflt/F");
    _hitdiag->Branch("active",&_active,"active/B");
    _hitdiag->Branch("use",&_use,"use/I");
    _hitdiag->Branch("edep",&_edep,"edep/F");
    _hitdiag->Branch("pdist",&_pdist,"pdist/F");
    _hitdiag->Branch("pperp",&_pperp,"pperp/F");
    _hitdiag->Branch("pmom",&_pmom,"pmom/F");
    _hitdiag->Branch("nmcsteps",&_nmcsteps,"nmcsteps/i");
    _hitdiag->Branch("mchpos",&_mchpos,"x/F:y/F:z/F");
    _hitdiag->Branch("mcrdrift",&_mcrdrift,"mcrdrift/F");
    _hitdiag->Branch("mchitt0",&_mchitt0,"mchitt0/F");
    _hitdiag->Branch("mcdmid",&_mcdmid,"mcdmid/F");

//      _hitdiag->Branch("shdir","CLHEP::HepVector",&_shdir);
    return _hitdiag;
  }

// find the MC truth objects in the event and set the local cache
  bool
  KalFitMC::findMCData(const art::Event& evt) {
    _mcdata.clear();
    _mchitsums.clear();
  // Get the persistent data about pointers to StepPointMCs
    art::Handle<PtrStepPointMCVectorCollection> mchitptrHandle;
    if(evt.getByLabel(_mcptrlabel,"StrawHitMCPtr",mchitptrHandle))
      _mcdata._mchitptr = mchitptrHandle.product();
  // Get the persistent data about the StepPointMCs, from the tracker and the virtual detectors
    art::Handle<StepPointMCCollection> mctrackerstepsHandle;
    if(evt.getByLabel(_mcstepslabel,"tracker",mctrackerstepsHandle))
      _mcdata._mcsteps = mctrackerstepsHandle.product();
    art::Handle<StepPointMCCollection> mcVDstepsHandle;
    if(evt.getByLabel(_mcstepslabel,"virtualdetector",mcVDstepsHandle))
      _mcdata._mcvdsteps = mcVDstepsHandle.product();
    art::Handle<SimParticleCollection> simParticlesHandle;
    if(evt.getByLabel(_simpartslabel,simParticlesHandle))
      _mcdata._simparts = simParticlesHandle.product();
    art::Handle<StrawDigiMCCollection> mcdigisHandle;
    if(evt.getByLabel(_mcdigislabel,"StrawHitMC",mcdigisHandle))
      _mcdata._mcdigis = mcdigisHandle.product();
// fill hit summary
    fillMCHitSummary();
    if( _mcdata.good()) {
// count # of conversion straw hits
      _strawhits = 0;
      _nchits = 0;
      _ncgood = 0;
      art::Handle<mu2e::StrawHitCollection> strawhitsH;
      if(evt.getByLabel(_strawhitslabel,strawhitsH))
	_strawhits = strawhitsH.product();
      if(_strawhits != 0){
	unsigned nstrs = _strawhits->size();
	for(unsigned istr=0; istr<nstrs;++istr){
	  const std::vector<MCHitSum>& mcsum = mcHitSummary(istr);
	  if(mcsum.size()>0){
	    bool conversion = (mcsum[0]._pdgid == 11 && mcsum[0]._gid == 2);
	    if(conversion){
	      ++_nchits;
	      if(mcsum[0]._spp->genParticle()->momentum().vect().mag()-mcsum[0]._mom.mag() < 5.0)
		++_ncgood;
	    }
	  }
	}
      }
      return true;
    }
    return false;
  }
  std::vector<int>const& KalFitMC::VDids(TRACKERPOS tpos) const {
    switch(tpos) {
      case trackerEnt: default:
	return _entvids;
      case trackerMid:
	return _midvids;
      case trackerExit:
	return _xitvids;
    }
  }

  double KalFitMC::MCT0(TRACKERPOS tpos) const {
    switch(tpos) {
      case trackerEnt: default:
	return _mcentinfo._time;
      case trackerMid:
	return _mcmidinfo._time;
      case trackerExit:
	return _mcxitinfo._time;
    }
  }
  double KalFitMC::MCMom(TRACKERPOS tpos) const {
    switch(tpos) {
      case trackerEnt: default:
	return _mcentinfo._mom;
      case trackerMid:
	return _mcmidinfo._mom;
      case trackerExit:
	return _mcxitinfo._mom;
    }
  }
  const helixpar& KalFitMC::MCHelix(TRACKERPOS tpos) const {
    switch(tpos) {
      case trackerEnt: default:
	return _mcentinfo._hpar;
      case trackerMid:
	return _mcmidinfo._hpar;
      case trackerExit:
	return _mcxitinfo._hpar;
    }
  }

  void
  KalFitMC::fillMCTrkInfo(MCStepItr const& imcs, MCTrkInfo& einfo) const {
    GlobalConstantsHandle<ParticleDataTable> pdt;
    GeomHandle<DetectorSystem> det;

    einfo._time = imcs->time();
    einfo._pdgid = imcs->simParticle()->pdgId(); 
    double charge = pdt->particle(imcs->simParticle()->pdgId()).ref().charge();
    CLHEP::Hep3Vector mom = imcs->momentum();
    einfo._mom = mom.mag();
    // need to transform into the tracker coordinate system
    CLHEP::Hep3Vector pos = det->toDetector(imcs->position());
    HepPoint ppos(pos.x(),pos.y(),pos.z());
    einfo._pos = pos;
    double hflt(0.0);
    HepVector parvec(5,0);
    GeomHandle<BFieldManager> bfmgr;
    CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(Hep3Vector(0.0,0.0,0.0));
    double bz = bfmgr->getBField(vpoint_mu2e).z();
    TrkHelixUtils::helixFromMom( parvec, hflt,ppos,
	mom,charge,bz);
    einfo._hpar = helixpar(parvec);
  }

  void
  KalFitMC::fillMCTrkInfo(SimParticle const& sp, MCTrkInfo& einfo) const {
    GlobalConstantsHandle<ParticleDataTable> pdt;
    GeomHandle<DetectorSystem> det;
    einfo._time = sp.startGlobalTime();
    einfo._pdgid = sp.pdgId();
    double charge = pdt->particle(sp.pdgId()).ref().charge();
    CLHEP::Hep3Vector mom = sp.startMomentum();
    einfo._mom = mom.mag();
    // need to transform into the tracker coordinate system
    CLHEP::Hep3Vector pos = det->toDetector(sp.startPosition());
    HepPoint ppos =(pos.x(),pos.y(),pos.z());
    einfo._pos = pos;
    double hflt(0.0);
    HepVector parvec(5,0);
    GeomHandle<BFieldManager> bfmgr;
    CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(Hep3Vector(0.0,0.0,0.0));
    double bz = bfmgr->getBField(vpoint_mu2e).z();
    TrkHelixUtils::helixFromMom( parvec, hflt,ppos,
	mom,charge,bz);
    einfo._hpar = helixpar(parvec);
  }
  
  KalFitMC::relation KalFitMC::relationship(StrawDigiMC const& mcd1, StrawDigiMC const& mcd2) {
    art::Ptr<SimParticle> ptr1, ptr2;
    if(mcd1.stepPointMC(StrawDigi::zero).isNonnull())
      ptr1 = mcd1.stepPointMC(StrawDigi::zero)->simParticle();
    else if(mcd1.stepPointMC(StrawDigi::one).isNonnull())
      ptr1 = mcd1.stepPointMC(StrawDigi::one)->simParticle();
    if(mcd2.stepPointMC(StrawDigi::zero).isNonnull())
      ptr2 = mcd2.stepPointMC(StrawDigi::zero)->simParticle();
    else if(mcd2.stepPointMC(StrawDigi::one).isNonnull())
      ptr2 = mcd2.stepPointMC(StrawDigi::one)->simParticle();
    return relationship(ptr1,ptr2);
  }

  KalFitMC::relation KalFitMC::relationship(art::Ptr<SimParticle> const& sppi,art::Ptr<SimParticle> const& sppj) {
    if(sppi.isNull() || sppj.isNull()) return none;
    if(sppi == sppj)return same;
    art::Ptr<SimParticle> pi = sppi->parent();
    art::Ptr<SimParticle> pj = sppj->parent();
    if(pi.isNonnull() && pi == sppj)return daughter;
    if(pj.isNonnull() && pj == sppi)return mother;
    if(pi.isNonnull() && pj.isNonnull()){
      if( pi == pj)return sibling;
      std::vector<art::Ptr<SimParticle> > pvi, pvj;
      pvi.push_back(sppi);
      pvj.push_back(sppj);
      while(pi.isNonnull()){
	pvi.push_back(pi);
	pi = pi->parent();
      }
      while(pj.isNonnull()){
	pvj.push_back(pj);
	pj = pj->parent();
      }
      std::vector<art::Ptr<SimParticle> >::iterator ifnd;
      ifnd = std::find(pvi.begin(),pvi.end(),sppj);
      if(ifnd != pvi.end())return udaughter;
      ifnd = std::find(pvj.begin(),pvj.end(),sppi);
      if(ifnd != pvj.end())return umother;
      for(size_t ii=0;ii<pvj.size();++ii){
	ifnd = std::find(pvi.begin(),pvi.end(),pvj[ii]);
	if(ifnd != pvi.end())return usibling;
      }
      for(size_t ii=0;ii<pvi.size();++ii){
	ifnd = std::find(pvj.begin(),pvj.end(),pvi[ii]);
	if(ifnd != pvj.end())return usibling;
      }
    }
    return none;
  }

}
