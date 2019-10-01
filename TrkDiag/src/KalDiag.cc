//
// MC functions associated with KalFit
//
#include "TrkDiag/inc/KalDiag.hh"
#include "TrkDiag/inc/TrkStrawHitInfo.hh"
#include "TrkDiag/inc/TrkCaloHitInfo.hh"
#include "TrkDiag/inc/TrkStrawHitInfoMC.hh"
//geometry
#include "GeometryService/inc/GeometryService.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "BFieldGeom/inc/BFieldConfig.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
// services
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "art_root_io/TFileService.h"
// data
#include "art/Framework/Principal/Event.h"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "RecoDataProducts/inc/TrkStrawHitSeed.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
// Utilities
// tracker
#include "TrackerGeom/inc/Straw.hh"
#include "TrkReco/inc/TrkUtilities.hh"
// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/TrkBase/TrkHelixUtils.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/TrkBase/TrkMomCalculator.hh"
#include "BTrk/BField/BFieldFixed.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/KalmanTrack/KalMaterial.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "Mu2eBTrk/inc/DetStrawElem.hh"
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
using CLHEP::Hep3Vector;

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
    _mcstepslabel(pset.get<string>("MCStepsLabel")),
    _mcstepsinstance(pset.get<string>("MCStepsInstance")),
    _simpartslabel(pset.get<string>("SimParticleLabel")),
    _simpartsinstance(pset.get<string>("SimParticleInstance")),
    _mcdigislabel(pset.get<string>("StrawDigiMCLabel")),
    _toff(pset.get<fhicl::ParameterSet>("TimeOffsets")),
    _fillmc(pset.get<bool>("FillMCInfo",true)),
    _debug(pset.get<int>("debugLevel",0)),
    _diag(pset.get<int>("diagLevel",1)),
    _uresid(pset.get<bool>("UnbiasedResiduals",true)),
    _mingood(pset.get<double>("MinimumGoodMomentumFraction",0.9)),
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
    // clear the branches
    reset();
    if(krep != 0){
  // fill track info for the tracker entrance
      fillTrkInfo(krep,_trkinfo);
    // compute hit information
      if(_diag > 1){
	fillHitInfo(krep,_tshinfo);
	fillMatInfo(krep,_tminfo);
      }
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
      } else if(_mcdata._simparts != 0) {
      // find 1st primary particle
	for ( auto isp = _mcdata._simparts->begin(); isp != _mcdata._simparts->end(); ++isp ){
	  if(isp->second.isPrimary()){
	    spp = art::Ptr<SimParticle>(_mcdata._simparthandle,isp->second.id().asInt());
	    break;
	  }
	}
      }
// Fill general MC information
      if(spp.isNonnull()){
	fillTrkInfoMC(spp,krep,_mcinfo);
// fill hit-specific MC information
	if(_diag > 1 && krep != 0)fillHitInfoMC(spp,krep,_tshinfomc);
// fill mc info at the particle production point
	fillTrkInfoMCStep(spp,_mcgeninfo);
      // find MC info at tracker
	if(_mcdata._mcvdsteps != 0){
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
      }
    }
    // fill the TTree if requested
    if(fill){
      _trkdiag->Fill();
    }
  }

  void
    KalDiag::fillTrkInfo(const KalRep* krep,TrkInfo& trkinfo) const {
    GeomHandle<VirtualDetector> vdg;
    GeomHandle<DetectorSystem> det;
//    if(krep != 0 && krep->fitCurrent()){
    if(krep != 0 ){
      trkinfo._status = krep->fitStatus().success();
      trkinfo._pdg = krep->particleType().particleType();
      trkinfo._t0 = krep->t0().t0();
      trkinfo._t0err = krep->t0().t0Err();
      trkinfo._nhits = krep->hitVector().size();
      trkinfo._ndof = krep->nDof();
      trkinfo._nactive = krep->nActive();
      trkinfo._chisq = krep->chisq();
      trkinfo._fitcon = TrkUtilities::chisqConsistency(krep);
      trkinfo._radlen = krep->radiationFraction();
      trkinfo._startvalid = krep->startValidRange();
      trkinfo._endvalid = krep->endValidRange();
// site counting
      int nmat(0), nmatactive(0), nbend(0);
      for(auto isite : krep->siteList()){
	if(isite->kalMaterial() != 0){
	  ++nmat;
	  if(isite->isActive())++nmatactive;
	}
	if(isite->kalBend() != 0)
	  ++nbend;
      }
      trkinfo._nmat = nmat;
      trkinfo._nmatactive = nmatactive;
      trkinfo._nbend = nbend;

//      Hep3Vector seedmom = TrkMomCalculator::vecMom(*(krep->seed()),krep->kalContext().bField(),0.0);
// trkinfo._seedmom = seedmom.mag();
      // count hits
      countHits(krep,trkinfo);
      // get the fit at the entrance to the tracker
      Hep3Vector entpos = det->toDetector(vdg->getGlobal(VirtualDetectorId::TT_FrontPA));
      double zent = entpos.z();
      // we don't know which way the fit is going: try both, and pick the one with the smallest flightlength
      double firsthitfltlen = krep->lowFitRange();
      double lasthitfltlen = krep->hiFitRange();
      double entlen = min(firsthitfltlen,lasthitfltlen);
      TrkHelixUtils::findZFltlen(krep->traj(),zent,entlen,0.1);
      // compute the tracker entrance fit information
      //      fillTrkFitInfo(krep,entlen,trkinfo._ent); // no longer stored in TrkInfo
    } else {
      // failed fit
      trkinfo._status = -krep->fitStatus().failure();
    }
  }

  void KalDiag::countHits(const KalRep* krep,TrkInfo& tinfo) const {
    // count number of hits with other (active) hits in the same panel
    tinfo._ndouble = tinfo._ndactive = tinfo._nnullambig = 0;
    tinfo._firstflt = tinfo._lastflt = 0.0;
    // loop over hits, and count
    TrkStrawHitVector tshv;
    convert(krep->hitVector(),tshv);
    bool first(false);
    tinfo._nhits = tshv.size();
    tinfo._nactive = tinfo._ndouble = tinfo._ndactive = tinfo._nnullambig = 0;
    for(auto ihit=tshv.begin(); ihit != tshv.end(); ++ihit) {
      const TrkStrawHit* tsh = *ihit;
      if(tsh != 0 ){
	if(tsh->isActive() && !first){
	  first = true;
	  tinfo._firstflt = tsh->fltLen();
	}
	if(tsh->isActive())tinfo._nactive++;
	if(tsh->ambig()==0)tinfo._nnullambig++;
	auto jhit = ihit; jhit++;
	if(jhit != tshv.end() && tsh->straw().id().uniquePanel() ==
	  (*jhit)->straw().id().uniquePanel()){
	    tinfo._ndouble++;
	    if(tsh->isActive())tinfo._ndactive++;
	}
      }
    }
    for(auto ihit = tshv.rbegin(); ihit != tshv.rend(); ++ihit){
      if((*ihit) != 0 && (*ihit)->isActive()){
	tinfo._lastflt = (*ihit)->fltLen();
	break;
      }
    }

 }

  void KalDiag::fillTrkFitInfo(const KalRep* krep, TrkFitInfo& trkfitinfo) const {
    GeomHandle<VirtualDetector> vdg;
    GeomHandle<DetectorSystem> det;
    // get the fit at the entrance to the tracker
    Hep3Vector entpos = det->toDetector(vdg->getGlobal(VirtualDetectorId::TT_FrontPA));
    double zent = entpos.z();
    // we don't know which way the fit is going: try both, and pick the one with the smallest flightlength
    double firsthitfltlen = krep->lowFitRange();
    double lasthitfltlen = krep->hiFitRange();
    double entlen = min(firsthitfltlen,lasthitfltlen);
    TrkHelixUtils::findZFltlen(krep->traj(),zent,entlen,0.1);
    // compute the tracker entrance fit information
    fillTrkFitInfo(krep,entlen,trkfitinfo);
  }

  void KalDiag::fillTrkFitInfo(const KalRep* krep,double fltlen,TrkFitInfo& trkfitinfo) const {
    //    trkfitinfo._fltlen = fltlen;
    // find momentum and parameters
    double loclen(0.0);
    const TrkSimpTraj* ltraj = krep->localTrajectory(fltlen,loclen);
    trkfitinfo._fitpar = helixpar(ltraj->parameters()->parameter());
    trkfitinfo._fitparerr = helixpar(ltraj->parameters()->covariance());
    Hep3Vector fitmom = krep->momentum(fltlen);
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
    TrkStrawHitVector tshv;
    convert(krep->hitVector(),tshv);
    for(auto ihit=tshv.begin(); ihit != tshv.end(); ++ihit) {
	const TrkStrawHit* tsh = *ihit;
	// loop over the hits and find the associated steppoints
	if(tsh != 0 && tsh->isActive()){
	  StrawDigiMC const& mcdigi = _mcdata._mcdigis->at(tsh->index());
	  art::Ptr<SimParticle> spp = mcdigi.earlyStepPointMC()->simParticle();
// see if this particle has already been found; if so, increment, if not, add it
	  bool found(false);
	  for(size_t isp=0;isp<sct.size();++isp){
// count direct daughter/parent as part the same particle
	    if(sct[isp]._spp == spp ){
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
    TrkStrawHitVector tshv;
    convert(krep->hitVector(),tshv);
    for(auto ihit=tshv.begin(); ihit != tshv.end(); ++ihit) {
      const TrkStrawHit* tsh = *ihit;
      if(tsh != 0){
        TrkStrawHitInfo tshinfo;
	fillHitInfo(tsh,tshinfo);
// count correlations with other TSH
	for(auto jhit=tshv.begin(); jhit != ihit; ++jhit){
	  const TrkStrawHit* otsh = *jhit;
	  if(otsh != 0){
	    if(tshinfo._plane ==  otsh->straw().id().getPlane() &&
		tshinfo._panel == otsh->straw().id().getPanel() ){
	      tshinfo._dhit = true;
	      if(otsh->isActive()){
		tshinfo._dactive = true;
		break;
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
    tshinfo._plane = tsh->straw().id().plane();
    tshinfo._panel = tsh->straw().id().panel();
    tshinfo._layer = tsh->straw().id().layer();
    tshinfo._straw = tsh->straw().id().straw();
    tshinfo._edep = tsh->comboHit().energyDep();
    // kludge CLHEP problem
    HepPoint hpos = tsh->hitTraj()->position(tsh->hitLen());
    tshinfo._poca = XYZVec(hpos.x(),hpos.y(),hpos.z());
    double resid,residerr;
    if(tsh->resid(resid,residerr,_uresid)){
      tshinfo._resid = resid;
      tshinfo._residerr = residerr;
    } else {
      tshinfo._resid = tshinfo._residerr = -1000.;
    }
    tshinfo._rdrift = tsh->driftRadius();
    tshinfo._rdrifterr = tsh->driftRadiusErr();
    tshinfo._wdist = tsh->comboHit().wireDist();
    tshinfo._werr = tsh->comboHit().wireRes();
    double rstraw = tsh->straw().getRadius();
    tshinfo._dx = sqrt(max(0.0,rstraw*rstraw-tshinfo._rdrift*tshinfo._rdrift));
    tshinfo._trklen = tsh->fltLen();
    tshinfo._hlen = tsh->hitLen();
    Hep3Vector tdir = tsh->trkTraj()->direction(tshinfo._trklen);
    tshinfo._wdot = tdir.dot(tsh->straw().getDirection());
    // for now approximate the local bfield direction as the z axis FIXME!!
    //    tshinfo._bdot = tdir.z();
    tshinfo._t0 = tsh->hitT0()._t0;
    // include signal propagation time correction
    tshinfo._ht = tsh->time()-tsh->signalTime();
//    tshinfo._dt = tsh->comboHit().dt();
    tshinfo._ambig = tsh->ambig();
    tshinfo._driftend = tsh->comboHit().driftEnd();
    tshinfo._tdrift = tsh->driftTime();
//    tshinfo._totcal = tsh->comboHit().TOT(StrawEnd::cal);
//    tshinfo._tothv = tsh->comboHit().TOT(StrawEnd::hv);
    if(tsh->hasResidual())
      tshinfo._doca = tsh->poca().doca();
    else
      tshinfo._doca = -100.0;
    //    tshinfo._exerr = tsh->driftVelocity()*tsh->temperature();
    //    tshinfo._penerr = tsh->penaltyErr();
    tshinfo._t0err = tsh->t0Err()/tsh->driftVelocity();
// cannot count correlations with other hits in this function; set to false
    tshinfo._dhit = tshinfo._dactive = false;
  }

  void KalDiag::fillCaloHitInfo(const TrkCaloHit* tch, TrkCaloHitInfo& tchinfo) const {
    tchinfo._active = tch->isActive();  
    tchinfo._did = tch->caloCluster().diskId();
    tchinfo._trklen = tch->fltLen();
    tchinfo._clen = tch->hitLen();
    HepPoint hpos = tch->hitTraj()->position(tch->hitLen());
    tchinfo._poca = XYZVec(hpos.x(),hpos.y(),hpos.z());
    if(tch->hasResidual())
      tchinfo._doca = tch->poca().doca();
    else
      tchinfo._doca = -100.0;
    tchinfo._t0 = tch->hitT0().t0();
    tchinfo._t0err = tch->hitT0().t0Err();
    tchinfo._ct = tch->time() - tch->timeOffset();
    tchinfo._cterr = tch->timeErr();
    tchinfo._edep = tch->caloCluster().energyDep();
  }

  void KalDiag::fillHitInfoMC(art::Ptr<SimParticle> const& pspp, const KalRep* krep,
     std::vector<TrkStrawHitInfoMC>& tshinfomcs) const {
    tshinfomcs.clear();
 // loop over hits
    TrkStrawHitVector tshv;
    convert(krep->hitVector(),tshv);
    for(auto ihit=tshv.begin(); ihit != tshv.end(); ++ihit) {
      const TrkStrawHit* tsh = *ihit;
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
    art::Ptr<StepPointMC> const& spmcp = mcdigi.earlyStepPointMC();
    art::Ptr<SimParticle> const& spp = spmcp->simParticle();
    // create MC info and fill
    tshinfomc._t0 = _toff.timeWithOffsetsApplied(*spmcp);
    tshinfomc._pdg = spp->pdgId();
    tshinfomc._proc = spp->originParticle().creationCode();
    tshinfomc._edep = mcdigi.energySum();
    tshinfomc._gen = -1;
    if(spp->genParticle().isNonnull())
      tshinfomc._gen = spp->genParticle()->generatorId().id();
    MCRelationship rel(pspp,spp);
    tshinfomc._rel = rel.relationship();
    // find the step midpoint
    Hep3Vector mcsep = spmcp->position()-straw.getMidPoint();
    Hep3Vector dir = spmcp->momentum().unit();
    tshinfomc._mom = spmcp->momentum().mag();
    tshinfomc._cpos = mcdigi.clusterPosition(mcdigi.earlyEnd());
    tshinfomc._len = mcsep.dot(straw.getDirection());
    Hep3Vector mcperp = (dir.cross(straw.getDirection())).unit();
    double dperp = mcperp.dot(mcsep);
    tshinfomc._twdot = dir.dot(straw.getDirection());
    tshinfomc._dist = fabs(dperp);
    tshinfomc._ambig = dperp > 0 ? -1 : 1; // follow TrkPoca convention
    // use 2-line POCA here
    TwoLinePCA pca(spmcp->position(),dir,straw.getMidPoint(),straw.getDirection());
    tshinfomc._doca = pca.dca();
  }

  void KalDiag::fillTrkInfoMC(art::Ptr<SimParticle> const&  spp,const KalRep* krep,
    TrkInfoMC& mcinfo) {
    // basic information
    if(spp->genParticle().isNonnull())
      mcinfo._gen = spp->genParticle()->generatorId().id();
    mcinfo._pdg = spp->pdgId();
    mcinfo._proc = spp->originParticle().creationCode();
    /*    art::Ptr<SimParticle> pp = spp->originParticle().parent();
    if(pp.isNonnull()){
      mcinfo._ppdg = pp->pdgId();
      mcinfo._pproc = pp->originParticle().creationCode();
      mcinfo._pmom = pp->startMomentum().vect().mag();
      if(pp->genParticle().isNonnull())
	mcinfo._pgen = pp->genParticle()->generatorId().id();
	}*/

    art::Ptr<SimParticle> sp = spp;
    // find the first parent which comes from a generator
    while(sp->genParticle().isNull() && sp->parent().isNonnull()){
      sp = sp->parent();
    }

    Hep3Vector mcmomvec = spp->startMomentum();
    double mcmom = mcmomvec.mag();
    // fill track-specific  MC info
    mcinfo._nactive = mcinfo._nhits = mcinfo._nambig = 0;
    if(krep != 0){
      TrkStrawHitVector tshv;
      convert(krep->hitVector(),tshv);
      for(auto ihit=tshv.begin(); ihit != tshv.end(); ++ihit) {
	const TrkStrawHit* tsh = *ihit;
	if(tsh != 0){
	  StrawDigiMC const& mcdigi = _mcdata._mcdigis->at(tsh->index());
	  art::Ptr<StepPointMC> const& spmcp = mcdigi.earlyStepPointMC();
	  if(spp == spmcp->simParticle()){
	    ++mcinfo._nhits;
	    // easiest way to get MC ambiguity is through info object
	    TrkStrawHitInfoMC tshinfomc;
	    fillHitInfoMC(spp,mcdigi,tsh->straw(),tshinfomc);
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
      if( imcd->stepPointMC(StrawEnd::cal)->simParticle() == spp){
	mcinfo._ndigi++;
	if(imcd->stepPointMC(StrawEnd::cal)->momentum().mag()/mcmom > _mingood)
	  mcinfo._ndigigood++;
      }
    }
  }

  void KalDiag::fillMatInfo(const KalRep* krep, std::vector<TrkStrawMatInfo>& tminfos ) const {
    tminfos.clear();
 // loop over sites, pick out the materials
    for(auto isite : krep->siteList()) {
      if(isite->kalMaterial() != 0){
	TrkStrawMatInfo tminfo;
	if(fillMatInfo(isite->kalMaterial(),tminfo))
	  tminfos.push_back(tminfo);
      }
    }
  }

  bool KalDiag::fillMatInfo(const KalMaterial* kmat, TrkStrawMatInfo& tminfo) const {
    bool retval(false);
    const DetStrawElem* delem = dynamic_cast<const DetStrawElem*>(kmat->detElem());
    if(delem != 0){
      retval = true;
      // KalMaterial info
      tminfo._active = kmat->isActive();
      tminfo._dp = kmat->momFraction();
      tminfo._radlen = kmat->radiationFraction();
      //      tminfo._sigMS = kmat->deflectRMS();
      // DetIntersection info
      const DetIntersection& dinter = kmat->detIntersection();
      //      tminfo._thit = (dinter.thit != 0);
      //      tminfo._thita = (dinter.thit != 0 && dinter.thit->isActive());
      tminfo._doca = dinter.dist;
      tminfo._tlen = dinter.pathlen;
      // straw information
      Straw const* straw = delem->straw();
      tminfo._plane = straw->id().getPlane();
      tminfo._panel = straw->id().getPanel();
      tminfo._layer = straw->id().getLayer();
      tminfo._straw = straw->id().getStraw();
    }
    return retval;
  }

  TTree*
  KalDiag::createTrkDiag(){
    art::ServiceHandle<art::TFileService> tfs;
    _trkdiag=tfs->make<TTree>("trkdiag","trk diagnostics");
    addBranches(_trkdiag);
    return _trkdiag;
  }

  void KalDiag::addBranches(TTree* trkdiag,const char* branchprefix){
    string brapre(branchprefix);
      // general track info
      trkdiag->Branch((brapre+"fit").c_str(),&_trkinfo,TrkInfo::leafnames().c_str());
// basic MC info
    if(_fillmc){
      // general MC info
      trkdiag->Branch((brapre+"mc").c_str(),&_mcinfo,TrkInfoMC::leafnames().c_str());
      // mc info at generation and several spots in the tracker
      trkdiag->Branch((brapre+"mcgen").c_str(),&_mcgeninfo,TrkInfoMCStep::leafnames().c_str());
      trkdiag->Branch((brapre+"mcent").c_str(),&_mcentinfo,TrkInfoMCStep::leafnames().c_str());
      trkdiag->Branch((brapre+"mcmid").c_str(),&_mcmidinfo,TrkInfoMCStep::leafnames().c_str());
      trkdiag->Branch((brapre+"mcxit").c_str(),&_mcxitinfo,TrkInfoMCStep::leafnames().c_str());
    }
// track hit info
    if(_diag > 1){
      trkdiag->Branch((brapre+"tsh").c_str(),&_tshinfo);
      if(_fillmc)trkdiag->Branch((brapre+"tshmc").c_str(),&_tshinfomc);
      trkdiag->Branch((brapre+"tm").c_str(),&_tminfo);
    }
  }

  // find the MC truth objects in the event and set the local cache
  bool
  KalDiag::findMCData(const art::Event& evt) {
    _mcdata.clear();
    if(_fillmc){
      // Get the persistent data about the StepPointMCs, from the tracker and the virtual detectors
      art::Handle<StepPointMCCollection> mcVDstepsHandle;
      if(evt.getByLabel(_mcstepslabel,_mcstepsinstance,mcVDstepsHandle))
	_mcdata._mcvdsteps = mcVDstepsHandle.product();
      if(evt.getByLabel(_simpartslabel,_simpartsinstance,_mcdata._simparthandle))
	_mcdata._simparts = _mcdata._simparthandle.product();
      art::Handle<StrawDigiMCCollection> mcdigisHandle;
      if(evt.getByLabel(_mcdigislabel,mcdigisHandle))
	_mcdata._mcdigis = mcdigisHandle.product();
      // update time offsets
      _toff.updateMap(evt);
      static bool first(true);
      if (first && !_mcdata.good()){
	first = false;
	std::cout <<"Warning: MC data is not complete as follows, proceeding" << std::endl;
	_mcdata.printPointerValues();
      }
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
    Hep3Vector mom = imcs->momentum();
    // need to transform into the tracker coordinate system
    Hep3Vector pos = det->toDetector(imcs->position());
    fillTrkInfoMCStep(mom,pos,charge,mcstepinfo);
  }

  void
  KalDiag::fillTrkInfoMCStep(art::Ptr<SimParticle> const& spp, TrkInfoMCStep& mcstepinfo) const {
    GlobalConstantsHandle<ParticleDataTable> pdt;
    GeomHandle<DetectorSystem> det;
    mcstepinfo._time = _toff.totalTimeOffset(spp) + spp->startGlobalTime();
    double charge = pdt->particle(spp->pdgId()).ref().charge();
    Hep3Vector mom = spp->startMomentum();
    // need to transform into the tracker coordinate system
    Hep3Vector pos = det->toDetector(spp->startPosition());
    fillTrkInfoMCStep(mom,pos,charge,mcstepinfo);
  }

  void
  KalDiag::fillTrkInfoMCStep(Hep3Vector const& mom, Hep3Vector const& pos, double charge, TrkInfoMCStep& mcstepinfo) const {
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;

    mcstepinfo._mom = Geom::toXYZVec(mom);
    mcstepinfo._pos = Geom::toXYZVec(pos);
    double hflt(0.0);
    HepVector parvec(5,0);
    static Hep3Vector vpoint_mu2e = det->toMu2e(Hep3Vector(0.0,0.0,0.0));
    static double bz = bfmgr->getBField(vpoint_mu2e).z();
    HepPoint ppos(pos.x(),pos.y(),pos.z());
    TrkHelixUtils::helixFromMom( parvec, hflt,ppos, mom,charge,bz);
    mcstepinfo._hpar = helixpar(parvec);
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

  unsigned KalDiag::countCEHits() const {
    unsigned ncehits = 0;
    unsigned nstrs = mcData()._mcdigis->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      StrawDigiMC const& mcdigi = mcData()._mcdigis->at(istr);
      art::Ptr<StepPointMC> const& spmcp = mcdigi.earlyStepPointMC();
      art::Ptr<SimParticle> const& spp = spmcp->simParticle();
      bool conversion = (spp->genParticle().isNonnull() && spp->genParticle()->generatorId().isConversion() && spmcp->momentum().mag()>90.0);
      if(conversion){
	++ncehits;
      }
    }
    return ncehits;
  }

}
