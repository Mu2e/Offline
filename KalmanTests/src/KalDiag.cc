//
// MC functions associated with KalFit
// $Id: KalDiag.cc,v 1.6 2014/09/22 12:13:17 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/09/22 12:13:17 $
//
//geometry
#include "GeometryService/inc/GeometryService.hh"
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
#include "KalmanTests/inc/KalDiag.hh"
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
#include <functional>

using namespace std; 
 
namespace mu2e 
{
  // comparison functor for ordering step points
  struct timecomp : public binary_function<MCStepItr,MCStepItr, bool> {
    bool operator()(MCStepItr x,MCStepItr y) { return x->time() < y->time(); }
  };

  struct spcountcomp : public binary_function<spcount, spcount , bool> {
    bool operator() (spcount a, spcount b) { return a._count > b._count; }
  };


  TrkStrawHitInfo::TrkStrawHitInfo() : _active(false), _usable(false), _device(-1),
  _sector(-1), _layer(-1), _straw(-1), _nplane(0), _npanel(0), _nlayer(0), 
  _z(-1000.0), _phi(-1000.0), _rho(-1000.0),
  _resid(-1000.0), _residerr(-1000.0), _rdrift(-1000.0), _rdrifterr(-1000.0),
  _trklen(-1000.0),_doca(-1000.0), _exerr(-1000.0), _penerr(-1000.0),
  _t0(-1000.0), _t0err(-1000.0), _ht(-1000.0), _tddist(-1000.0), _tdderr(-1000.0),
  _hlen(-1000.0), _edep(-1000.0), _dx(-1000.0), _ambig(-1)
  {}

  TrkStrawHitInfoMC::TrkStrawHitInfoMC() : 
  _mcpdg(-1), _mcgen(-1), _mcproc(-1),
  _mct0(-1000.0), _mcht(-1000.0), _mcdist(-1000.0), _mclen(-1000.0),
  _mcedep(-1000.0),_mcr(-1000.0),_mcphi(-1000.0),
  _mcambig(-100), _xtalk(false) {}

  KalDiag::~KalDiag(){}
  
  KalDiag::KalDiag(fhicl::ParameterSet const& pset) :
    _mcptrlabel(pset.get<string>("MCPtrLabel")),
    _mcstepslabel(pset.get<string>("MCStepsLabel")),
    _simpartslabel(pset.get<string>("SimParticleLabel")),
    _simpartsinstance(pset.get<string>("SimParticleInstance")),
    _mcdigislabel(pset.get<string>("StrawHitMCLabel")),
    _toff(pset.get<fhicl::ParameterSet>("TimeOffsets")),
    _fillmc(pset.get<bool>("FillMCInfo",true)),
    _debug(pset.get<int>("debugLevel",0)),
    _diag(pset.get<int>("diagLevel",1)),
    _uresid(pset.get<bool>("UnbiasedResiduals",true)),
    _mingood(pset.get<double>("MinimumGoodMomentumFraction",0.9)),
    _mintrkmom(pset.get<double>("minTrkMom",60.0)),
    _mct0err(pset.get<double>("mcT0Err",0.1)),
    _mcambig(pset.get<bool>("mcAmbiguity",false)),
    _minnhits(pset.get<unsigned>("minNHits",10)),
    _maxnhits(pset.get<unsigned>("maxNHits",120)),
    _purehits(pset.get<bool>("pureHits",false)),
    _trkdiag(0)
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
  KalDiag::findMCSteps(StepPointMCCollection const* mcsteps, cet::map_vector_key const& trkid, vector<int> const& vids,
  vector<MCStepItr>& steps) {
    steps.clear();
    if(mcsteps != 0){
      // Loop over the step points, and find the one corresponding to the given detector
      for( MCStepItr imcs =mcsteps->begin();imcs!= mcsteps->end();imcs++){
	if(vids.size() == 0 ||  (imcs->trackId() == trkid && find(vids.begin(),vids.end(),imcs->volumeId()) != vids.end())){
	  steps.push_back(imcs);
	}
      }
      // sort these in time
      sort(steps.begin(),steps.end(),timecomp());
    }
  }

  void
  KalDiag::kalDiag(const KalRep* krep,bool fill) {
    GeomHandle<VirtualDetector> vdg;
    GeomHandle<DetectorSystem> det;
    if(_trkdiag == 0)createTrkDiag();
    // clear the branches
    reset();
   // null pointer to SimParticle
    art::Ptr<SimParticle> spp;;
    if(krep != 0) {
      if(krep->fitCurrent())
	_fitstatus = krep->fitStatus().success();
      else
    // failed fit
	_fitstatus = -krep->fitStatus().failure();
      _t0 = krep->t0().t0();
      _t0err = krep->t0().t0Err();
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
      double entlen = min(firsthitfltlen,lasthitfltlen);
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
// add MC info
      if(_fillmc){
	vector<spcount> sct;
	findMCTrk(krep,sct);
	_nmc = sct.size();
	if(sct.size()>0 && sct[0]._spp.isNonnull()){
	  spp = sct[0]._spp;
	  fillMCTrkInfo(spp);
	}
      }
      // compute hit information no matter what, as summary information is still used      if(_diag > 1)
      hitsDiag(krep,spp);
    } else if(_fillmc && _mcdata._simparts != 0){
// Assume the 1st particle is the primary 
      for ( auto isp = _mcdata._simparts->begin(); isp != _mcdata._simparts->end(); ++isp ){
	if(isp->second.isPrimary()){
	  art::Ptr<SimParticle> spp(_mcdata._simparthandle,isp->second.id().asInt());
	  fillMCTrkInfo(spp);
	  break;
	}
      }
    }
  // fill the TTree if requested
    if(fill)_trkdiag->Fill(); 
  }

  void KalDiag::findMCTrk(const KalRep* krep,art::Ptr<SimParticle>& spp) {
    static art::Ptr<SimParticle> nullvec;
    spp = nullvec;
    vector<spcount> sct;
    findMCTrk(krep,sct);
    if(sct.size()>0)
      spp = sct[0]._spp;
  }

  
  void
  KalDiag::findMCTrk(const KalRep* krep,vector<spcount>& sct) {
    sct.clear();
// find the SimParticles which contributed hits.
    if(_mcdata._mcdigis != 0) {
// get the straw hits from the track
      const TrkHotList* hots = krep->hotList();
      for(TrkHotList::hot_iterator ihot=hots->begin();ihot != hots->end();++ihot){
	const TrkStrawHit* tsh = dynamic_cast<const TrkStrawHit*>(ihot.get());
	// loop over the hits and find the associated steppoints
	if(tsh != 0 && tsh->isActive()){
	  StrawDigiMC const& mcdigi = _mcdata._mcdigis->at(tsh->index());
	  StrawDigi::TDCChannel itdc = StrawDigi::zero;
	  if(!mcdigi.hasTDC(StrawDigi::zero))
	    itdc = StrawDigi::one;
	  art::Ptr<SimParticle> spp = mcdigi.stepPointMC(itdc)->simParticle();
// see if this particle has already been found; if so, increment, if not, add it
	  bool found(false);
	  for(size_t isp=0;isp<sct.size();++isp){
// count direct daughter/parent as part the same particle
	    if(sct[isp]._spp == spp ){
//	    || sct[isp]._spp->parent()==spp || spp->parent() == sct[isp]._spp){
	      found = true;
	      sct[isp].append(spp);
	      break;
	    }
	  }
	  if(!found)sct.push_back(spp);
	}
      }
    }
    // sort by # of contributions
    sort(sct.begin(),sct.end(),spcountcomp());
  }
  
  void KalDiag::hitsDiag(const KalRep* krep,art::Ptr<SimParticle> const& pspp) {
    _tshinfo.clear();
    _tshinfomc.clear();
    _nmcactive = _nmchits = _nmcgood = 0;
    _nmcambig = 0;
    _ndouble = _ndactive = 0;
    bool mcgood = _mcdata.good();
    double mcmom(100.0);
    if(pspp.isNonnull())
      mcmom = pspp->startMomentum().mag();
 // loop over hits
    const TrkHotList* hots = krep->hotList();
    for(TrkHotList::hot_iterator ihot=hots->begin();ihot != hots->end();++ihot){
      const TrkStrawHit* tsh = dynamic_cast<const TrkStrawHit*>(ihot.get());
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
	if(tsh->resid(resid,residerr,_uresid)){
	  tshinfo._resid = resid;
	  tshinfo._residerr = residerr;
	} else {
	  tshinfo._resid = tshinfo._residerr = -1000.;
	}
	tshinfo._rdrift = tsh->driftRadius();
	tshinfo._rdrifterr = tsh->driftRadiusErr();
	double rstraw = tsh->straw().getRadius();
	tshinfo._dx = sqrt(max(0.0,rstraw*rstraw-tshinfo._rdrift*tshinfo._rdrift));
	tshinfo._trklen = tsh->fltLen();
	tshinfo._hlen = tsh->hitLen();
	Hep3Vector tdir = tsh->trkTraj()->direction(tshinfo._trklen);
	tshinfo._wdot = tdir.dot(tsh->straw().getDirection());
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
// count correlations with other TSH
	tshinfo._nlayer = tshinfo._npanel = tshinfo._nplane = 0;
	bool dhit(false);
	bool dactive(false);
	for(TrkHotList::hot_iterator jhot=hots->begin();jhot != hots->end();++jhot){
	  if(jhot != ihot){
	    const TrkStrawHit* otsh = dynamic_cast<const TrkStrawHit*>(jhot.get());
	    if(otsh != 0){
	      if(tshinfo._device ==  otsh->straw().id().getDevice()){
		tshinfo._nplane++;
		if(tshinfo._sector == otsh->straw().id().getSector()){
		  tshinfo._npanel++;
		  dhit = true;
		  if(tsh->isActive() && otsh->isActive())dactive = true;
		  if(tshinfo._layer == otsh->straw().id().getLayer()){
		    tshinfo._nlayer++;
		  }
		}
	      }
	    }
	  }
	}
	if(dhit)_ndouble++;
	if(dactive)_ndactive++;
	_tshinfo.push_back(tshinfo);
	// MC information
	if(_fillmc && mcgood){
	  StrawDigiMC const& mcdigi = _mcdata._mcdigis->at(tsh->index());
	  // use TDC channel 0 to define the MC match
	  StrawDigi::TDCChannel itdc = StrawDigi::zero;
	  if(!mcdigi.hasTDC(StrawDigi::one)) itdc = StrawDigi::one;
	  art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
	  art::Ptr<SimParticle> const& spp = spmcp->simParticle();
	  // create MC info and fill
	  TrkStrawHitInfoMC tshinfomc;
	  tshinfomc._mct0 = _toff.timeWithOffsetsApplied(*spmcp);
	  tshinfomc._mcht = mcdigi.wireEndTime(itdc);
	  tshinfomc._mcpdg = spp->pdgId();
	  tshinfomc._mcproc = spp->creationCode();
	  tshinfomc._mcedep = mcdigi.energySum();
	  tshinfomc._mcgen = -1;
	  if(spp->genParticle().isNonnull())
	    tshinfomc._mcgen = spp->genParticle()->generatorId().id();
	  tshinfomc._mcrel = relationship(pspp,spp);
	  // find the step midpoint
	  Hep3Vector mcsep = spmcp->position()-tsh->straw().getMidPoint();
	  Hep3Vector dir = spmcp->momentum().unit();
	  Hep3Vector mcperp = (dir.cross(tsh->straw().getDirection())).unit();
	  double dperp = mcperp.dot(mcsep);
	  tshinfomc._mcr =spmcp->position().perp();
	  tshinfomc._mcphi =spmcp->position().phi();
	  tshinfomc._mcdist = fabs(dperp);
	  tshinfomc._mcambig = dperp > 0 ? -1 : 1; // follow TrkPoca convention
	  tshinfomc._mclen = mcsep.dot(tsh->straw().getDirection());
	  tshinfomc._xtalk = spmcp->strawIndex() != tsh->strawHit().strawIndex();
	  // count hits from the primary particle
	  if(spp == pspp){
	    ++_nmchits;
	    if(tsh->isActive())++_nmcactive;
	    // count hits with at least givien fraction of the original momentum as 'good'
	    if(spmcp->momentum().mag()/mcmom > _mingood )++_nmcgood;
	    // count hits with correct left-right ambiguity
	    if(tshinfo._active && tshinfomc._mcambig*tshinfo._ambig > 0)++_nmcambig;
	  }
// fill entry
	  _tshinfomc.push_back(tshinfomc);
	}
      }
    }
  }

  void KalDiag::fillMCTrkInfo(art::Ptr<SimParticle> const&  spp) {
    GeomHandle<VirtualDetector> vdg;
    GeomHandle<DetectorSystem> det;
// basic information
    if(spp->genParticle().isNonnull())
      _mcgenid = spp->genParticle()->generatorId().id();
    _mcpdgid = spp->pdgId();
    _mcproc = spp->creationCode();
    art::Ptr<SimParticle> pp = spp->parent();
    if(pp.isNonnull()){
      _mcppdgid = pp->pdgId();
      _mcpproc = pp->creationCode();
      if(pp->genParticle().isNonnull())
	_mcpgenid = pp->genParticle()->generatorId().id();
    }

// find the mc info at the entrance to the detector
    fillMCTrkInfo(spp,_mcinfo);
    // find MC info at tracker
    cet::map_vector_key trkid = spp->id();
    vector<MCStepItr> entsteps,midsteps,xitsteps;
    findMCSteps(_mcdata._mcvdsteps,trkid,_entvids,entsteps);
    if(entsteps.size() > 0 && vdg->exist(entsteps[0]->volumeId()))
      fillMCTrkInfo(entsteps.front(),_mcentinfo);
    findMCSteps(_mcdata._mcvdsteps,trkid,_midvids,midsteps);
    if(midsteps.size() > 0 && vdg->exist(midsteps[0]->volumeId()))
      fillMCTrkInfo(midsteps.front(),_mcmidinfo);
    findMCSteps(_mcdata._mcvdsteps,trkid,_xitvids,xitsteps);
    if(xitsteps.size() > 0 && vdg->exist(xitsteps[0]->volumeId()))
      fillMCTrkInfo(xitsteps.front(),_mcxitinfo);

// count the # of tracker hits (digis) generated by this particle
    _npdigi = _npdgood = 0;
    for(auto imcd = _mcdata._mcdigis->begin(); imcd !=_mcdata._mcdigis->end();++imcd){
      if( imcd->hasTDC(StrawDigi::zero) && imcd->stepPointMC(StrawDigi::zero)->simParticle() == spp){
      _npdigi++;
      if(imcd->stepPointMC(StrawDigi::zero)->momentum().mag()/spp->startMomentum().mag() > _mingood)
	_npdgood++;
      }
    }
  }

// map daughters onto parents within the context of an associated set of StepPointMCs (like from a StrawHit).
  void KalDiag::findRelatives(PtrStepPointMCVector const& mcptr,map<SPPtr,SPPtr>& mdmap){
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
          map<SPPtr,SPPtr>::iterator ifnd = mdmap.find(sppp);
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
  KalDiag::createTrkDiag(){
    art::ServiceHandle<art::TFileService> tfs;
    _trkdiag=tfs->make<TTree>("trkdiag","trk diagnostics");
    _trkdiag->Branch("fitstatus",&_fitstatus,"fitstatus/I");
    _trkdiag->Branch("t0",&_t0,"t0/F");
    _trkdiag->Branch("t0err",&_t0err,"t0err/F");
    _trkdiag->Branch("nhits",&_nhits,"nhits/I");
    _trkdiag->Branch("ndof",&_ndof,"ndof/I");
    _trkdiag->Branch("niter",&_niter,"niter/I");
    _trkdiag->Branch("nactive",&_nactive,"nactive/I");
    _trkdiag->Branch("ndouble",&_ndouble,"ndouble/I");
    _trkdiag->Branch("ndactive",&_ndactive,"ndactive/I");
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
    if(_fillmc){
      _trkdiag->Branch("npdigi",&_npdigi,"npdigi/I");
      _trkdiag->Branch("npdgood",&_npdgood,"npdgood/I");
      _trkdiag->Branch("nmc",&_nmc,"nmc/I");
      _trkdiag->Branch("nmcactive",&_nmcactive,"nmcactive/I");
      _trkdiag->Branch("nmchits",&_nmchits,"nmchits/I");
      _trkdiag->Branch("nmcgood",&_nmcgood,"nmcgood/I");
      _trkdiag->Branch("nmcambig",&_nmcambig,"nmcambig/I");
      _trkdiag->Branch("mcpdgid",&_mcpdgid,"mcpdgid/I");
      _trkdiag->Branch("mcgenid",&_mcgenid,"mcgenid/I");
      _trkdiag->Branch("mcproc",&_mcproc,"mcproc/I");
      _trkdiag->Branch("mcppdgid",&_mcppdgid,"mcppdgid/I");
      _trkdiag->Branch("mcpgenid",&_mcpgenid,"mcpgenid/I");
      _trkdiag->Branch("mcpproc",&_mcpproc,"mcpproc/I");
      // mc info at production and several spots in the tracker
      _trkdiag->Branch("mcinfo",&_mcinfo,"mct0/F:mcmom/F:mcx/F:mcy/F:mcz/F:mcd0/F:mcp0/F:mcom/F:mcz0/F:mctd/F");
      _trkdiag->Branch("mcentinfo",&_mcentinfo,"mcentt0/F:mcentmom/F:mcentx/F:mcenty/F:mcentz/F:mcentd0/F:mcentp0/F:mcentom/F:mcentz0/F:mcenttd/F");
      _trkdiag->Branch("mcmidinfo",&_mcmidinfo,"mcmidt0/F:mcmidmom/F:mcmidx/F:mcmidy/F:mcmidz/F:mcmidd0/F:mcmidp0/F:mcmidom/F:mcmidz0/F:mcmidtd/F");
      _trkdiag->Branch("mcxitinfo",&_mcxitinfo,"mcxitt0/F:mcxitmom/F:mcxitx/F:mcxity/F:mcxitz/F:mcxitd0/F:mcxitp0/F:mcxitom/F:mcxitz0/F:mcxittd/F");
    }
// track hit info    
    if(_diag > 1){
      _trkdiag->Branch("tshinfo",&_tshinfo);
      if(_fillmc)_trkdiag->Branch("tshinfomc",&_tshinfomc);
    }
    return _trkdiag;
  }

  // find the MC truth objects in the event and set the local cache
  bool
  KalDiag::findMCData(const art::Event& evt) {
    _mcdata.clear();
    if(_fillmc){
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
      if(evt.getByLabel(_simpartslabel,_simpartsinstance,_mcdata._simparthandle))
	_mcdata._simparts = _mcdata._simparthandle.product();
      art::Handle<StrawDigiMCCollection> mcdigisHandle;
      if(evt.getByLabel(_mcdigislabel,"StrawHitMC",mcdigisHandle))
	_mcdata._mcdigis = mcdigisHandle.product();
      // update time offsets
      _toff.updateMap(evt);
      return _mcdata.good();
    }
    return true;
  }

  vector<int>const& KalDiag::VDids(TRACKERPOS tpos) const {
    switch(tpos) {
      case trackerEnt: default:
	return _entvids;
      case trackerMid:
	return _midvids;
      case trackerExit:
	return _xitvids;
    }
  }

  const helixpar& KalDiag::MCHelix(TRACKERPOS tpos) const {
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
  KalDiag::fillMCTrkInfo(MCStepItr const& imcs, MCTrkInfo& einfo) const {
    GlobalConstantsHandle<ParticleDataTable> pdt;
    GeomHandle<DetectorSystem> det;
    einfo._time = _toff.timeWithOffsetsApplied(*imcs);
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
    TrkHelixUtils::helixFromMom( parvec, hflt,ppos, mom,charge,bz);
    einfo._hpar = helixpar(parvec);
  }

  void
  KalDiag::fillMCTrkInfo(art::Ptr<SimParticle> const& spp, MCTrkInfo& einfo) const {
    GlobalConstantsHandle<ParticleDataTable> pdt;
    GeomHandle<DetectorSystem> det;
//    einfo._time = spp->startGlobalTime();
//    TimeOffset interface forces a copy
    art::Ptr<SimParticle> sp(spp);
    einfo._time = _toff.totalTimeOffset(sp) + spp->startGlobalTime();
    double charge = pdt->particle(spp->pdgId()).ref().charge();
    CLHEP::Hep3Vector mom = spp->startMomentum();
    einfo._mom = mom.mag();
    // need to transform into the tracker coordinate system
    CLHEP::Hep3Vector pos = det->toDetector(spp->startPosition());
    HepPoint ppos =(pos.x(),pos.y(),pos.z());
    einfo._pos = pos;
    double hflt(0.0);
    HepVector parvec(5,0);
    GeomHandle<BFieldManager> bfmgr;
    CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(Hep3Vector(0.0,0.0,0.0));
    double bz = bfmgr->getBField(vpoint_mu2e).z();
    TrkHelixUtils::helixFromMom( parvec, hflt,ppos, mom,charge,bz);
    einfo._hpar = helixpar(parvec);
  }
 
 // StrawDigiMCs record the StepPoint that pushed the electronics over threshold; use that
 // define the relationship

  KalDiag::relation KalDiag::relationship(StrawDigiMC const& mcd1, StrawDigiMC const& mcd2) {
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

  KalDiag::relation KalDiag::relationship(art::Ptr<SimParticle> const& sppi,art::Ptr<SimParticle> const& sppj) {
    if(sppi.isNull() || sppj.isNull()) return none;
    if(sppi == sppj)return same;
    art::Ptr<SimParticle> pi = sppi->parent();
    art::Ptr<SimParticle> pj = sppj->parent();
    if(pi.isNonnull() && pi == sppj)return daughter;
    if(pj.isNonnull() && pj == sppi)return mother;
    if(pi.isNonnull() && pj.isNonnull()){
      if( pi == pj)return sibling;
      vector<art::Ptr<SimParticle> > pvi, pvj;
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
      vector<art::Ptr<SimParticle> >::iterator ifnd;
      ifnd = find(pvi.begin(),pvi.end(),sppj);
      if(ifnd != pvi.end())return udaughter;
      ifnd = find(pvj.begin(),pvj.end(),sppi);
      if(ifnd != pvj.end())return umother;
      for(size_t ii=0;ii<pvj.size();++ii){
	ifnd = find(pvi.begin(),pvi.end(),pvj[ii]);
	if(ifnd != pvi.end())return usibling;
      }
      for(size_t ii=0;ii<pvi.size();++ii){
	ifnd = find(pvj.begin(),pvj.end(),pvi[ii]);
	if(ifnd != pvj.end())return usibling;
      }
    }
    return none;
  }

  void KalDiag::reset() {
    // reset ttree variables that might otherwise be stale and misleading
    _fitstatus = -1000;
    _nhits = _nactive = _ndouble = _ndactive = _ndof = _niter = _nsites = -1;
    _fitmom = _seedmom = _t0 = _t0err = -1.0;
    _chisq = _fitcon = _radlen = _firstflt = _lastflt = -1.0;
    _tshinfo.clear();
    if(_fillmc){
      _nmc = _nmcactive = _nmchits = _nmcgood = _nmcambig = -1;
      _mcpdgid = _mcgenid = _mcproc = -1;
      _mcppdgid = _mcpgenid = _mcpproc = -1;
      _mcinfo.reset();
      _mcentinfo.reset();
      _mcmidinfo.reset();
      _mcxitinfo.reset();
      _tshinfomc.clear();
    }
  }

 // define seed helix, t0, and hits coming from a given particle using MC truth.  Note that the input
// trkdef object must reference a valid straw hit collection
  bool
  KalDiag::trkFromMC(cet::map_vector_key const& trkid,TrkDef& mytrk) {
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
  

}
