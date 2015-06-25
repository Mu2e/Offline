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
#include "BTrk/BaBar/BaBar.hh"
#include "KalmanTests/inc/KalDiag.hh"
#include "BTrk/TrkBase/TrkHelixUtils.hh"
#include "BTrk/TrkBase/TrkHotList.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/TrkBase/TrkMomCalculator.hh"
#include "BTrk/BField/BFieldFixed.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
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
#include "TString.h"
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
// initialize TrkQual MVA.  Note the weight file is passed in from the KalDiag config
    fhicl::ParameterSet mvapset = pset.get<fhicl::ParameterSet>("TrkQualMVA",fhicl::ParameterSet());
    mvapset.put<string>("MVAWeights",pset.get<string>("TrkQualWeights","KalmanTests/test/TrkQual.weights.xml"));
    _trkqualmva.reset(new MVATools(mvapset));
    _trkqualmva->initMVA();
    if(_debug>0)_trkqualmva->showMVA();
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
    // clear the branches
    reset();
    if(krep != 0){
  // fill track info for the tracker entrance
      fillTrkInfo(krep,_trkinfo);
    // compute hit information 
      if(_diag > 1)fillHitInfo(krep,_tshinfo);
    }
// if requested, add MC info
    if(_fillmc){
   // null pointer to SimParticle
      art::Ptr<SimParticle> spp;;
      if(krep != 0){
  // if the fit exists, find the MC particle which best matches it
	vector<spcount> sct;
	findMCTrk(krep,sct);
	if(sct.size()>0 && sct[0]._spp.isNonnull()){
	  spp = sct[0]._spp;
	}
      } else { 
      // find 1st primary particle
	for ( auto isp = _mcdata._simparts->begin(); isp != _mcdata._simparts->end(); ++isp ){
	  if(isp->second.isPrimary()){
	    spp = art::Ptr<SimParticle>(_mcdata._simparthandle,isp->second.id().asInt());
	    break;
	  }
	}
      }
// Fill general MC information
      fillTrkInfoMC(spp,krep,_mcinfo);
// fill hit-specific MC information
      if(_diag > 1 && krep != 0)fillHitInfoMC(spp,krep,_tshinfomc);
// fill mc info at the particle production point
      fillTrkInfoMCStep(spp,_mcgeninfo);
      // find MC info at tracker
      cet::map_vector_key trkid = spp->id();
      vector<MCStepItr> entsteps,midsteps,xitsteps;
      findMCSteps(_mcdata._mcvdsteps,trkid,_entvids,entsteps);
      if(entsteps.size() > 0 && vdg->exist(entsteps[0]->volumeId()))
	fillTrkInfoMCStep(entsteps.front(),_mcentinfo);
      findMCSteps(_mcdata._mcvdsteps,trkid,_midvids,midsteps);
      if(midsteps.size() > 0 && vdg->exist(midsteps[0]->volumeId()))
	fillTrkInfoMCStep(midsteps.front(),_mcmidinfo);
      findMCSteps(_mcdata._mcvdsteps,trkid,_xitvids,xitsteps);
      if(xitsteps.size() > 0 && vdg->exist(xitsteps[0]->volumeId()))
	fillTrkInfoMCStep(xitsteps.front(),_mcxitinfo);
    }
    // fill the TTree if requested
    if(fill){
      _trkdiag->Fill();
    }
  }

  void
  KalDiag::fillTrkInfo(const KalRep* krep,TrkInfo& trkinfo) {
    GeomHandle<VirtualDetector> vdg;
    GeomHandle<DetectorSystem> det;
    if(krep->fitCurrent())
      trkinfo._fitstatus = krep->fitStatus().success();
    else
      // failed fit
      trkinfo._fitstatus = -krep->fitStatus().failure();
    trkinfo._fitpart = krep->particleType().particleType();
    trkinfo._t0 = krep->t0().t0();
    trkinfo._t0err = krep->t0().t0Err();
    trkinfo._nhits = krep->hotList()->nHit();
    trkinfo._ndof = krep->nDof();
    trkinfo._nactive = krep->nActive();
    trkinfo._chisq = krep->chisq();
    trkinfo._fitcon = krep->chisqConsistency().significanceLevel();
    trkinfo._radlen = krep->radiationFraction();
    trkinfo._firstflt = krep->firstHit()->globalLength();
    trkinfo._lastflt = krep->lastHit()->globalLength();
    CLHEP::Hep3Vector seedmom = TrkMomCalculator::vecMom(*(krep->seed()),krep->kalContext().bField(),0.0);
    trkinfo._seedmom = seedmom.mag();
// count # of double hits
    countDoubles(krep,trkinfo._ndouble, trkinfo._ndactive);
   // get the fit at the entrance to the tracker
    CLHEP::Hep3Vector entpos = det->toDetector(vdg->getGlobal(VirtualDetectorId::TT_FrontPA));
    double zent = entpos.z();
    // we don't know which way the fit is going: try both, and pick the one with the smallest flightlength
    double firsthitfltlen = krep->lowFitRange(); 
    double lasthitfltlen = krep->hiFitRange();
    double entlen = min(firsthitfltlen,lasthitfltlen);
    TrkHelixUtils::findZFltlen(krep->traj(),zent,entlen,0.1);
    // compute the tracker entrance fit information
    fillTrkFitInfo(krep,entlen,trkinfo._ent);
// use the above information to compute the TrkQual value.
    fillTrkQual(trkinfo); 
  } 

  void
  KalDiag::countDoubles(const KalRep* krep,int& ndouble, int& ndactive) const {
// count number of hits with other (active) hits in the same panel
    ndouble = ndactive = 0;
// loop over hits, and count
    const TrkHotList* hots = krep->hotList();
    for(TrkHotList::hot_iterator ihot=hots->begin();ihot != hots->end();++ihot){
      const TrkStrawHit* tsh = dynamic_cast<const TrkStrawHit*>(ihot.get());
      if(tsh != 0){
	bool isdouble(false);
	bool dactive(false);
// count correlations with other TSH
	for(TrkHotList::hot_iterator jhot=hots->begin();jhot != hots->end();++jhot){
	  if(jhot != ihot){
	    const TrkStrawHit* otsh = dynamic_cast<const TrkStrawHit*>(jhot.get());
	    if(otsh != 0){
	      if(tsh->straw().id().getDevice() ==  otsh->straw().id().getDevice() &&
		  tsh->straw().id().getSector() == otsh->straw().id().getSector() ){
		  isdouble = true;
		if(otsh->isActive()){
		  dactive = true;
		  break;
		}
	      }
	    }
	  }
	}
	if(isdouble)ndouble++;
	if(dactive)ndactive++;
      }
    }
  }

  void
  KalDiag::fillTrkFitInfo(const KalRep* krep,double fltlen,TrkFitInfo& trkfitinfo) {
  // find momentum and parameters
    double loclen(0.0);
    const TrkSimpTraj* ltraj = krep->localTrajectory(fltlen,loclen);
    trkfitinfo._fitpar = helixpar(ltraj->parameters()->parameter());
    trkfitinfo._fitparerr = helixpar(ltraj->parameters()->covariance());
    CLHEP::Hep3Vector fitmom = krep->momentum(fltlen);
    BbrVectorErr momerr = krep->momentumErr(fltlen);
    trkfitinfo._fitmom = fitmom.mag();
    Hep3Vector momdir = fitmom.unit();
    HepVector momvec(3);
    for(int icor=0;icor<3;icor++)
      momvec[icor] = momdir[icor];
    trkfitinfo._fitmomerr = sqrt(momerr.covMatrix().similarity(momvec));
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
  
  void KalDiag::fillHitInfo(const KalRep* krep, std::vector<TrkStrawHitInfo>& tshinfos ) const { 
    tshinfos.clear();
 // loop over hits
    const TrkHotList* hots = krep->hotList();
    for(TrkHotList::hot_iterator ihot=hots->begin();ihot != hots->end();++ihot){
      const TrkStrawHit* tsh = dynamic_cast<const TrkStrawHit*>(ihot.get());
      if(tsh != 0){
        TrkStrawHitInfo tshinfo;
	fillHitInfo(tsh,tshinfo);
// count correlations with other TSH
	for(TrkHotList::hot_iterator jhot=hots->begin();jhot != hots->end();++jhot){
	  if(jhot != ihot){
	    const TrkStrawHit* otsh = dynamic_cast<const TrkStrawHit*>(jhot.get());
	    if(otsh != 0){
	      if(tshinfo._device ==  otsh->straw().id().getDevice() &&
		  tshinfo._sector == otsh->straw().id().getSector() ){
		tshinfo._dhit = true;
		if(otsh->isActive()){
		  tshinfo._dactive = true;
		  break;
		}
	      }
	    }
	  }
	}
	tshinfos.push_back(tshinfo);
      }
    }
  }

  void KalDiag::fillHitInfo(const TrkStrawHit* tsh,TrkStrawHitInfo& tshinfo) const {
    tshinfo._active = tsh->isActive();
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
// cannot count correlations with other hits in this function; set to false
    tshinfo._dhit = tshinfo._dactive = false;
  }

  void KalDiag::fillHitInfoMC(art::Ptr<SimParticle> const& pspp, const KalRep* krep, 
     std::vector<TrkStrawHitInfoMC>& tshinfomcs) const {
    tshinfomcs.clear();
 // loop over hits
    const TrkHotList* hots = krep->hotList();
    for(TrkHotList::hot_iterator ihot=hots->begin();ihot != hots->end();++ihot){
      const TrkStrawHit* tsh = dynamic_cast<const TrkStrawHit*>(ihot.get());
      if(tsh != 0){
	StrawDigiMC const& mcdigi = _mcdata._mcdigis->at(tsh->index());
	TrkStrawHitInfoMC tshinfomc;
	fillHitInfoMC(pspp,mcdigi,tsh->straw(),tshinfomc);
	tshinfomcs.push_back(tshinfomc);
      }
    }
  }

  void KalDiag::fillHitInfoMC(art::Ptr<SimParticle> const& pspp, StrawDigiMC const& mcdigi,Straw const& straw, 
    TrkStrawHitInfoMC& tshinfomc) const {
    // use TDC channel 0 to define the MC match
    StrawDigi::TDCChannel itdc = StrawDigi::zero;
    if(!mcdigi.hasTDC(itdc)) itdc = StrawDigi::one;
    art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
    art::Ptr<SimParticle> const& spp = spmcp->simParticle();
    // create MC info and fill
    tshinfomc._t0 = _toff.timeWithOffsetsApplied(*spmcp);
    tshinfomc._ht = mcdigi.wireEndTime(itdc);
    tshinfomc._pdg = spp->pdgId();
    tshinfomc._proc = spp->creationCode();
    tshinfomc._edep = mcdigi.energySum();
    tshinfomc._gen = -1;
    if(spp->genParticle().isNonnull())
      tshinfomc._gen = spp->genParticle()->generatorId().id();
    tshinfomc._rel = relationship(pspp,spp);
    // find the step midpoint
    Hep3Vector mcsep = spmcp->position()-straw.getMidPoint();
    Hep3Vector dir = spmcp->momentum().unit();
    tshinfomc._mom = spmcp->momentum().mag();
    tshinfomc._r =spmcp->position().perp();
    tshinfomc._phi =spmcp->position().phi();
    Hep3Vector mcperp = (dir.cross(straw.getDirection())).unit();
    double dperp = mcperp.dot(mcsep);
    tshinfomc._dist = fabs(dperp);
    tshinfomc._ambig = dperp > 0 ? -1 : 1; // follow TrkPoca convention
    tshinfomc._len = mcsep.dot(straw.getDirection());
    tshinfomc._xtalk = spmcp->strawIndex() != mcdigi.strawIndex();
  }

  void KalDiag::fillTrkInfoMC(art::Ptr<SimParticle> const&  spp,const KalRep* krep,
    TrkInfoMC& mcinfo) {
    // basic information
    if(spp->genParticle().isNonnull())
      mcinfo._gen = spp->genParticle()->generatorId().id();
    mcinfo._pdg = spp->pdgId();
    mcinfo._proc = spp->creationCode();
    art::Ptr<SimParticle> pp = spp->parent();
    if(pp.isNonnull()){
      mcinfo._ppdg = pp->pdgId();
      mcinfo._pproc = pp->creationCode();
      if(pp->genParticle().isNonnull())
	mcinfo._pgen = pp->genParticle()->generatorId().id();
    }
    CLHEP::Hep3Vector mcmomvec = spp->startMomentum();
    double mcmom = mcmomvec.mag(); 
    // fill track-specific  MC info
    mcinfo._nactive = mcinfo._nhits = mcinfo._ngood = mcinfo._nambig = 0;
    if(krep != 0){
      const TrkHotList* hots = krep->hotList();
      for(TrkHotList::hot_iterator ihot=hots->begin();ihot != hots->end();++ihot){
	const TrkStrawHit* tsh = dynamic_cast<const TrkStrawHit*>(ihot.get());
	if(tsh != 0){
	  StrawDigiMC const& mcdigi = _mcdata._mcdigis->at(tsh->index());
	  StrawDigi::TDCChannel itdc = StrawDigi::zero;
	  if(!mcdigi.hasTDC(StrawDigi::one)) itdc = StrawDigi::one;
	  art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
	  if(spp == spmcp->simParticle()){
	    ++mcinfo._nhits;
	    // easiest way to get MC ambiguity is through info object
	    TrkStrawHitInfoMC tshinfomc;
	    fillHitInfoMC(spp,mcdigi,tsh->straw(),tshinfomc);
	    // count hits with at least givien fraction of the original momentum as 'good'
	    if(tshinfomc._mom/mcmom > _mingood )++mcinfo._ngood;
	    if(tsh->isActive()){
	      ++mcinfo._nactive;
	    // count hits with correct left-right iguity
	      if(tsh->ambig()*tshinfomc._ambig > 0)++mcinfo._nambig;
	    }
	  }
	}
      }
    }

    // count the # of tracker hits (digis) generated by this particle
    mcinfo._ndigi = mcinfo._ndigigood = 0;
    for(auto imcd = _mcdata._mcdigis->begin(); imcd !=_mcdata._mcdigis->end();++imcd){
      if( imcd->hasTDC(StrawDigi::zero) && imcd->stepPointMC(StrawDigi::zero)->simParticle() == spp){
	mcinfo._ndigi++;
	if(imcd->stepPointMC(StrawDigi::zero)->momentum().mag()/spp->startMomentum().mag() > _mingood)
	  mcinfo._ndigigood++;
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
  KalDiag::createTrkDiag(TTree* trkdiag,const char* branchprefix){
    string brapre(branchprefix);
    if(trkdiag == 0){
      art::ServiceHandle<art::TFileService> tfs;
      _trkdiag=tfs->make<TTree>("trkdiag","trk diagnostics");
    } else 
      _trkdiag = trkdiag;
      // general track info
      _trkdiag->Branch((brapre+"fit").c_str(),&_trkinfo,TrkInfo::leafnames().c_str());
// basic MC info
    if(_fillmc){
      // general MC info
      _trkdiag->Branch((brapre+"mc").c_str(),&_mcinfo,TrkInfoMC::leafnames().c_str());
      // mc info at generation and several spots in the tracker
      _trkdiag->Branch((brapre+"mcgen").c_str(),&_mcgeninfo,TrkInfoMCStep::leafnames().c_str());
      _trkdiag->Branch((brapre+"mcent").c_str(),&_mcentinfo,TrkInfoMCStep::leafnames().c_str());
      _trkdiag->Branch((brapre+"mcmid").c_str(),&_mcmidinfo,TrkInfoMCStep::leafnames().c_str());
      _trkdiag->Branch((brapre+"mcxit").c_str(),&_mcxitinfo,TrkInfoMCStep::leafnames().c_str());
    }
// track hit info    
    if(_diag > 1){
      _trkdiag->Branch((brapre+"tsh").c_str(),&_tshinfo);
      if(_fillmc)_trkdiag->Branch((brapre+"tshmc").c_str(),&_tshinfomc);
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
      if (!_mcdata.good()) _mcdata.printPointerValues();
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
  KalDiag::fillTrkInfoMCStep(MCStepItr const& imcs, TrkInfoMCStep& mcstepinfo) const {
    GlobalConstantsHandle<ParticleDataTable> pdt;
    GeomHandle<DetectorSystem> det;
    mcstepinfo._time = _toff.timeWithOffsetsApplied(*imcs);
    double charge = pdt->particle(imcs->simParticle()->pdgId()).ref().charge();
    CLHEP::Hep3Vector mom = imcs->momentum();
    // need to transform into the tracker coordinate system
    CLHEP::Hep3Vector pos = det->toDetector(imcs->position());
    fillTrkInfoMCStep(mom,pos,charge,mcstepinfo);
  }

  void
  KalDiag::fillTrkInfoMCStep(art::Ptr<SimParticle> const& spp, TrkInfoMCStep& mcstepinfo) const {
    GlobalConstantsHandle<ParticleDataTable> pdt;
    GeomHandle<DetectorSystem> det;
    mcstepinfo._time = _toff.totalTimeOffset(spp) + spp->startGlobalTime();
    double charge = pdt->particle(spp->pdgId()).ref().charge();
    CLHEP::Hep3Vector mom = spp->startMomentum();
    // need to transform into the tracker coordinate system
    CLHEP::Hep3Vector pos = det->toDetector(spp->startPosition());
    fillTrkInfoMCStep(mom,pos,charge,mcstepinfo);
  }

  void
  KalDiag::fillTrkInfoMCStep(CLHEP::Hep3Vector const& mom, CLHEP::Hep3Vector const& pos, double charge, TrkInfoMCStep& mcstepinfo) const {
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;

    mcstepinfo._mom = mom.mag();
    mcstepinfo._pos = pos;
    double hflt(0.0);
    HepVector parvec(5,0);
    static CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(Hep3Vector(0.0,0.0,0.0));
    static double bz = bfmgr->getBField(vpoint_mu2e).z();
    HepPoint ppos(pos.x(),pos.y(),pos.z());
    TrkHelixUtils::helixFromMom( parvec, hflt,ppos, mom,charge,bz);
    mcstepinfo._hpar = helixpar(parvec);
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

  void KalDiag::fillTrkQual(TrkInfo& trkinfo) const {
    static std::vector<double> trkqualvec; // input variables for TrkQual computation
    trkqualvec.resize(8);
    trkqualvec[0] = trkinfo._nactive; // # of active hits
    trkqualvec[1] = (float)trkinfo._nactive/(float)trkinfo._nhits;  // Fraction of active hits
    trkqualvec[2] = log10(trkinfo._fitcon); // fit chisquared consistency
    trkqualvec[3] = trkinfo._ent._fitmomerr; // estimated momentum error
    trkqualvec[4] = trkinfo._t0err;  // estimated t0 error
    trkqualvec[5] = trkinfo._ent._fitpar._d0; // d0 value
    trkqualvec[6] = trkinfo._ent._fitpar._d0+2.0/trkinfo._ent._fitpar._om; // maximum radius of fit
    trkqualvec[7] = (float)trkinfo._ndactive/(float)trkinfo._nactive;  // fraction of double hits (2 or more in 1 panel)
    trkinfo._trkqual = _trkqualmva->evalMVA(trkqualvec);
  }
  
  void KalDiag::reset() {
    // reset ttree variables that might otherwise be stale and misleading
    _trkinfo.reset();
    _tshinfo.clear();
    if(_fillmc){
      _mcinfo.reset();
      _mcgeninfo.reset();
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
