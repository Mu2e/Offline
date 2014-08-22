//
// Look for particles coming from the calorimeter and reflecting back in the
// magnetic mirror
//
// $Id: Reflect_module.cc,v 1.11 2014/08/22 19:55:50 brownd Exp $
// $Author: brownd $
// $Date: 2014/08/22 19:55:50 $
//
// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
// services
#include "BFieldGeom/inc/BFieldConfig.hh"
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "MCDataProducts/inc/VirtualDetectorId.hh"
#include "GeometryService/inc/DetectorSystem.hh"
// data
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "RecoDataProducts/inc/TrkExtTrajCollection.hh"
// ROOT incldues
#include "TTree.h"
#include "TH1F.h"
// Need this for the BaBar headers.
using namespace CLHEP;
// BaBar includes
#include "BaBar/BaBar.hh"
#include "KalmanTrack/KalRep.hh"
#include "KalmanTrack/KalHit.hh"
#include "TrkBase/TrkParticle.hh"
#include "TrkBase/TrkHelixUtils.hh"
// mu2e tracking
#include "KalmanTests/inc/TrkFitDirection.hh"
#include "KalmanTests/inc/KalDiag.hh"
// C++ includes.
#include <iostream>
#include <string>
#include <math.h>
// This is fragile and needs to be last until CLHEP is
// properly qualified and included in the BaBar classes.
#include "KalmanTests/inc/KalRepCollection.hh"

using namespace std;

namespace mu2e {
// struct to describe a fit.  Use root types as this
// gets put into a TTree
  struct FitInfo {
    Int_t _fitpart; // particle type
    Int_t _fitstatus; // fit status
    Int_t _nactive; // # active hits
    Float_t _chisq,_fitcon; // fit chisq, consistency
    Float_t _t0, _t0err; // t0 and error (at tracker midplane)
    Float_t _entmom, _entmomerr; // momentum at tracker entrance
    helixpar _entpar, _entperr; // parameters and errors and the tracker entrance
    threevec _entpdir, _entpos; // entrance momentum direciton and position
    Float_t _tent, _tenterr; // time and error at entrance to tracker
    Float_t _entflt; // flight length at entrance
  };

  class Reflect : public art::EDAnalyzer {
  public:

    explicit Reflect(fhicl::ParameterSet const& pset);
    virtual ~Reflect() { }

    void beginJob();
    void analyze(const art::Event& e);

  private:

    // Module label of the module that performed the fits.
    std::string _fitterPrefix;
    std::string _extModName;
    std::vector<std::string> _udname, _umname; // upstream particle data product names
    std::vector<std::string> _ddname, _dmname; // downstream particle data product names
    int _ipart; // particle PDG code
    double _mindmom, _maxdmom, _maxdtd;
    double _maxdp0, _mindt0, _maxdt0;
    double _zent; // z position of the entrance of the tracker
    unsigned _eventid;
    bool _extrapolate;
    // diagnostic of Kalman fit
    KalDiag _kdiag;
    // TTree for studying reflecting fits
    TTree* _reflect;
    // TTree branches
    Int_t _upart, _dpart, _ufitstat, _dfitstat;
    Int_t _unactive, _dnactive;
    Float_t _uchisq, _dchisq, _ufitcon, _dfitcon;
    Float_t _ut0, _ut0err, _dt0, _dt0err;
    Float_t _umom, _umomerr, _dmom, _dmomerr;
    helixpar _upar, _uperr, _dpar, _dperr;
    threevec _upos, _dpos, _updir, _dpdir;
    Float_t _utent, _dtent, _utenterr, _dtenterr;
    Float_t _uentf, _dentf;
    MCTrkInfo _umcinfo,_dmcinfo;
    Float_t _et0;
    threevec _emom, _epos;
    Float_t _pt0;
    Int_t _ppdg;
    threevec _ppos, _pmom, _pppos;
    Int_t _uextnpa,_dextnpa,_uextnst,_dextnst;
    Float_t _uextdppa,_dextdppa,_uextdpst,_dextdpst;
// create 
    void createTree();
// fill tree
    void fillTree(FitInfo const& uinfo, FitInfo const& dinfo);
// fill fit information
    void fillFitInfo(const KalRep* krep,FitInfo& fitinfo) const;
    void fillParentInfo(art::Ptr<SimParticle> sp); 
    // Function to pair upstream and downstream fits
    bool reflection(FitInfo const& uinfo, FitInfo const& dinfo) const;
    void getEntranceZ();
  };

  Reflect::Reflect(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset),
    _fitterPrefix(pset.get<string>("fitterPrefix","TPR")),
    _extModName(pset.get<string>("ExtrapolationModuleName","TrkExt")),
    _ipart(pset.get<int>("ParticleCode",11)),
    _mindmom(pset.get<double>("MinDeltaMom",-10.0)),
    _maxdmom(pset.get<double>("MaxDeltaMom",10.0)),
    _maxdtd(pset.get<double>("MaxDeltaTanDip",0.25)),
    _maxdp0(pset.get<double>("MaxDeltaPhi0",0.25)),
    _mindt0(pset.get<double>("MinDeltaT0",-50)),
    _maxdt0(pset.get<double>("MaxDeltaT0",-110)),
    _extrapolate(pset.get<bool>("Extrapolate",true)),
    _kdiag(pset.get<fhicl::ParameterSet>("KalDiag",fhicl::ParameterSet())),
    _reflect(0)
  {
// construct the data product instance names for particles 
    TrkFitDirection udir(TrkFitDirection::upstream);
    TrkFitDirection ddir(TrkFitDirection::downstream);
    for(int isign=1;isign>-2;isign-=2){
      TrkParticle tpart((TrkParticle::type)(_ipart*isign));
      std::string udname = udir.name() + tpart.name();
      _udname.push_back(udname);
      std::string umname = _fitterPrefix + udname;
      _umname.push_back(umname);
      std::string ddname = ddir.name() + tpart.name();
      _ddname.push_back(ddname);
      std::string dmname = _fitterPrefix + ddname;
      _dmname.push_back(dmname);
    }
  }

  void Reflect::beginJob( ){
// create the tree
    createTree();
// initialize counter
    _eventid = 0;
  }

  void Reflect::getEntranceZ( ){
 // get the virtual detector at the tracker entrance and take its z position
    GeomHandle<VirtualDetector> vdg;
    GeomHandle<DetectorSystem> det;
    CLHEP::Hep3Vector entpos = det->toDetector(vdg->getGlobal(VirtualDetectorId::TT_FrontPA));
    _zent = entpos.z();
  } 

  // For each event, look at tracker hits and calorimeter hits.
  void Reflect::analyze(const art::Event& event) {
    if(_eventid==0)getEntranceZ();
    _eventid++;
// get MC info
    bool hasmc = _kdiag.findMCData(event);
// loop over particle type
    for(size_t ie=0;ie<_udname.size();++ie){
      TrkExtTrajCollection const* uext(0);
      TrkExtTrajCollection const* dext(0);
 
      art::Handle<KalRepCollection> utrksHandle;
      event.getByLabel(_umname[ie],_udname[ie],utrksHandle);
      if(!utrksHandle.isValid())continue;
      KalRepCollection const& utrks = *utrksHandle;
      // now downstream
      art::Handle<KalRepCollection> dtrksHandle;
      event.getByLabel(_dmname[ie],_ddname[ie],dtrksHandle);
      if(!dtrksHandle.isValid())continue;
      KalRepCollection const& dtrks = *dtrksHandle;
// if we're including extrapolation information, find it in the event
      if(_extrapolate){
	art::Handle<TrkExtTrajCollection> utrkextHandle;
	event.getByLabel(_extModName,_udname[ie],utrkextHandle);
	art::Handle<TrkExtTrajCollection> dtrkextHandle;
	event.getByLabel(_extModName,_ddname[ie],dtrkextHandle);
	// look for extrapolation information
	if(utrkextHandle.isValid() && dtrkextHandle.isValid()){
	  uext = &*utrkextHandle;
	  dext = &*dtrkextHandle;
	}
      }
// loop over pairs, and see if any match.
      for ( size_t iue=0; iue < utrks.size(); ++iue ){
	KalRep const* ukrep = utrks.get(iue);
	if ( ukrep != 0 && ukrep->fitCurrent() ){
	  FitInfo uinfo;
	  fillFitInfo(ukrep,uinfo);
	  for ( size_t ide=0; ide < dtrks.size(); ++ide ){
	    KalRep const* dkrep = dtrks.get(ide);
	    if ( dkrep != 0 && dkrep->fitCurrent() ){
	      FitInfo dinfo;
	      fillFitInfo(dkrep,dinfo);
	      if(reflection(uinfo,dinfo)){
// fill the branches with this info
		fillTree(uinfo,dinfo);
		_umcinfo = _dmcinfo = MCTrkInfo();
// get MC info for the upstream and downstream tracks
		if(hasmc){
		  art::Ptr<SimParticle> umcinfo, dmcinfo;
		  _kdiag.findMCTrk(ukrep,umcinfo);
		  _kdiag.findMCTrk(dkrep,dmcinfo);
// use these to find the points where the true particle enters the tracker
		  if(umcinfo.isNonnull() && dmcinfo.isNonnull() && 
		    umcinfo == dmcinfo){
		    std::vector<MCStepItr> steps;
		    _kdiag.findMCSteps(_kdiag.mcData()._mcvdsteps,umcinfo->id(),_kdiag.VDids(KalDiag::trackerEnt),steps);
		    if(steps.size() == 2){
// These are sorted by time: first should be upstream, second down
		      _kdiag.fillMCTrkInfo(steps[0],_umcinfo);
		      _kdiag.fillMCTrkInfo(steps[1],_dmcinfo);
		      if(_extrapolate && uext != 0 && dext != 0){
		        TrkExtTraj const& utrkext = (*uext)[iue];
		        TrkExtTraj const& dtrkext = (*dext)[ide];
			_uextnpa = utrkext.getNPAHits(); 
			_dextnpa = dtrkext.getNPAHits(); 
			_uextnst = utrkext.getNSTHits(); 
			_dextnst = dtrkext.getNSTHits(); 
			_uextdppa = utrkext.getDeltapPA(); 
			_dextdppa = dtrkext.getDeltapPA(); 
			_uextdpst = utrkext.getDeltapST(); 
			_dextdpst = dtrkext.getDeltapST(); 
		      } 
// fill info about the electron origin and the parent of this electron (IE the cosmic muon)
		      fillParentInfo(umcinfo);
		    } else
		      std::cout << "Didn't find 2 steps" << std::endl;
		  } else
		    std::cout << "MC info doesn't match " << std::endl;
		} else
		  std::cout << "No MC info " << std::endl;
		_reflect->Fill();
	      }
	    }
	  }
	}
      }
    }
  }

// fill fit information from a KalRep
  void
  Reflect::fillFitInfo(const KalRep* krep,FitInfo& fitinfo) const {
    fitinfo._fitpart = krep->particleType().particleType();
    fitinfo._fitstatus = krep->fitStatus().success();
    fitinfo._nactive = krep->nActive();
    fitinfo._chisq = krep->chisq();
    fitinfo._fitcon = krep->chisqConsistency().significanceLevel();
    fitinfo._t0 = krep->t0().t0();
    fitinfo._t0err = krep->t0().t0Err();
    // find where the track enters the tracker, and get its parameters
    // there
    double firsthitfltlen = krep->firstHit()->kalHit()->hitOnTrack()->fltLen() - 10;
    double lasthitfltlen = krep->lastHit()->kalHit()->hitOnTrack()->fltLen() - 10;
    double entfltlen = std::min(firsthitfltlen,lasthitfltlen);
    TrkHelixUtils::findZFltlen(krep->traj(),_zent,entfltlen,0.1); 
    double loclen(0.0);
    const TrkSimpTraj* ltraj = krep->localTrajectory(entfltlen,loclen);
    fitinfo._entpar = helixpar(ltraj->parameters()->parameter());
    fitinfo._entperr = helixpar(ltraj->parameters()->covariance());
    HepPoint epos = krep->position(entfltlen);
    fitinfo._entpos = Hep3Vector(epos.x(),epos.y(),epos.z());
    BbrVectorErr momerr = krep->momentumErr(entfltlen);
    CLHEP::Hep3Vector fitmom = krep->momentum(entfltlen);
    double entmom = fitmom.mag();
    fitinfo._entmom = entmom;
    Hep3Vector momdir = fitmom.unit();
    fitinfo._entpdir = momdir;
    HepVector momvec(3);
    for(int icor=0;icor<3;icor++)
      momvec[icor] = momdir[icor];
    fitinfo._entmomerr = sqrt(momerr.covMatrix().similarity(momvec));
    fitinfo._entflt = entfltlen;
// estimate the time at tracker entrance, assuming a constant average momentum
    double flt0(0.0);
    TrkHelixUtils::findZFltlen(krep->traj(),0.0,flt0,0.1);
    double dflt = entfltlen-flt0; 
    double beta = krep->particleType().beta(entmom);
    double tflt = dflt/(beta*CLHEP::c_light);
    double mom0 = krep->momentum(flt0).mag();
    double beta0 = krep->particleType().beta(mom0);
    double tflt0 = dflt/(beta0*CLHEP::c_light);
    fitinfo._tent = 0.5*(tflt+tflt0) + krep->t0().t0();
// estimate the time error as the t0 error plus the time due to momentum change
    double dt = tflt-tflt0;
    fitinfo._tenterr = sqrt(krep->t0().t0Err()*krep->t0().t0Err() + dt*dt/12.0);

  }
  bool 
  Reflect::reflection(FitInfo const& uinfo, FitInfo const& dinfo) const {
// check basic info
    bool retval = uinfo._fitstatus > 0 && dinfo._fitstatus > 0 &&
      uinfo._fitpart == dinfo._fitpart;
    if(retval){
      // rough cuts on momentum and direction at tracker entrance
      double dmom = uinfo._entmom-dinfo._entmom;
      double dtd = uinfo._entpar._td + dinfo._entpar._td; // note sign change!!!
      double dp0 = uinfo._entpar._p0 - dinfo._entpar._p0;
      if(fabs(dp0) > M_PI){
	if(dp0 > 0.0)
	  dp0 -= 2*M_PI;
	else
	  dp0 += 2*M_PI;
      }
      double dt = uinfo._t0 - dinfo._t0;
      retval &= dmom > _mindmom && dmom < _maxdmom && fabs(dtd) < _maxdtd &&
	fabs(dp0) < _maxdp0 && dt < _mindt0 && dt > _maxdt0;
    }
    return retval;
  }

  void
  Reflect::createTree() {
    art::ServiceHandle<art::TFileService> tfs;
    _reflect=tfs->make<TTree>("reflect","reflection diagnostics");
// create the branches
    _reflect->Branch("eventid",&_eventid,"eventid/I");
    // upstream information
    _reflect->Branch("upart",&_upart,"upart/I");
    _reflect->Branch("ufitstat",&_ufitstat,"ufitstat/I");
    _reflect->Branch("unactive",&_unactive,"unactive/I");
    _reflect->Branch("uchisq",&_uchisq,"uchisq/F");
    _reflect->Branch("ufitcon",&_ufitcon,"ufitcon/F");
    _reflect->Branch("ut0",&_ut0,"ut0/F");
    _reflect->Branch("ut0err",&_ut0err,"ut0err/F");
    _reflect->Branch("umom",&_umom,"umom/F");
    _reflect->Branch("umomerr",&_umomerr,"umomerr/F");
    _reflect->Branch("upar",&_upar,"ud0/F:up0/F:uom/F:uz0/F:utd/F");
    _reflect->Branch("uperr",&_uperr,"ud0err/F:up0err/F:uomerr/F:uz0err/F:utderr/F");
    _reflect->Branch("upos",&_upos,"ux/F:uy/F:uz/F");
    _reflect->Branch("updir",&_updir,"upx/F:upy/F:upz/F");
    _reflect->Branch("utent",&_utent,"utent/F");
    _reflect->Branch("utenterr",&_utenterr,"utenterr/F");
    _reflect->Branch("uentf",&_uentf,"uentf/F");
    _reflect->Branch("umcinfo",&_umcinfo,"umcpdgid/I:umctime/F:umcmom/F:umcx/F:umcy/F:umcz/F:umcd0/F:umcp0/F:umcom/F:umcz0/F:umctd/F");
    // downstream information
    _reflect->Branch("dpart",&_dpart,"dpart/I");
    _reflect->Branch("dfitstat",&_dfitstat,"dfitstat/I");
    _reflect->Branch("dnactive",&_dnactive,"dnactive/I");
    _reflect->Branch("dchisq",&_dchisq,"dchisq/F");
    _reflect->Branch("dfitcon",&_dfitcon,"dfitcon/F");
    _reflect->Branch("dt0",&_dt0,"dt0/F");
    _reflect->Branch("dt0err",&_dt0err,"dt0err/F");
    _reflect->Branch("dmom",&_dmom,"dmom/F");
    _reflect->Branch("dmomerr",&_dmomerr,"dmomerr/F");
    _reflect->Branch("dpar",&_dpar,"dd0/F:dp0/F:dom/F:dz0/F:dtd/F");
    _reflect->Branch("dperr",&_dperr,"dd0err/F:dp0err/F:domerr/F:dz0err/F:dtderr/F");
    _reflect->Branch("dpos",&_dpos,"dx/F:dy/F:dz/F");
    _reflect->Branch("dpdir",&_dpdir,"dpx/F:dpy/F:dpz/F");
    _reflect->Branch("dtent",&_dtent,"dtent/F");
    _reflect->Branch("dtent",&_dtent,"dtent/F");
    _reflect->Branch("dtenterr",&_dtenterr,"dtenterr/F");
    _reflect->Branch("dentf",&_dentf,"dentf/F");
    _reflect->Branch("dmcinfo",&_dmcinfo,"dmcpdgid/I:dmctime/F:dmcmom/F:dmcx/F:dmcy/F:dmcz/F:dmcd0/F:dmcp0/F:dmcom/F:dmcz0/F:dmctd/F");
    // general information about production electron and muon
    _reflect->Branch("emom",&_emom,"emx/F:emy/F:emz/F");
    _reflect->Branch("et0",&_et0,"et0/F");
    _reflect->Branch("epos",&_epos,"ex/F:ey/F:ez/F");
    _reflect->Branch("ppdg",&_ppdg,"ppdg/I");
    _reflect->Branch("pt0",&_pt0,"pt0/F");
    _reflect->Branch("ppos",&_ppos,"px/F:py/F:pz/F");
    _reflect->Branch("pmom",&_pmom,"pmx/F:pmy/F:pmz/F");
    _reflect->Branch("pppos",&_pppos,"ppx/F:ppy/F:ppz/F");
    if(_extrapolate){
      _reflect->Branch("uextnPA",&_uextnpa,"uextnpa/I");
      _reflect->Branch("dextnPA",&_dextnpa,"dextnpa/I");
      _reflect->Branch("uextnST",&_uextnst,"uextnst/I");
      _reflect->Branch("dextnST",&_dextnst,"dextnst/I");
      _reflect->Branch("uextdpPA",&_uextdppa,"uextdppa/F");
      _reflect->Branch("dextdpPA",&_dextdppa,"dextdppa/F");
      _reflect->Branch("uextdpST",&_uextdpst,"uextdpst/F");
      _reflect->Branch("dextdpST",&_dextdpst,"dextdpst/F");
    }
  }

  void
  Reflect::fillTree(FitInfo const& uinfo, FitInfo const& dinfo) {
    _upart = uinfo._fitpart;
    _ufitstat = uinfo._fitstatus;
    _unactive = uinfo._nactive;
    _uchisq = uinfo._chisq;
    _ufitcon = uinfo._fitcon;
    _ut0 = uinfo._t0;
    _ut0err = uinfo._t0err;
    _umom = uinfo._entmom;
    _umomerr = uinfo._entmomerr;
    _upar = uinfo._entpar;
    _uperr = uinfo._entperr;
    _upos = uinfo._entpos;
    _updir = uinfo._entpdir;
    _utent = uinfo._tent;
    _utenterr = uinfo._tenterr;
    _uentf = uinfo._entflt;
// downstream info
    _dpart = dinfo._fitpart;
    _dfitstat = dinfo._fitstatus;
    _dnactive = dinfo._nactive;
    _dchisq = dinfo._chisq;
    _dfitcon = dinfo._fitcon;
    _dt0 = dinfo._t0;
    _dt0err = dinfo._t0err;
    _dmom = dinfo._entmom;
    _dmomerr = dinfo._entmomerr;
    _dpar = dinfo._entpar;
    _dperr = dinfo._entperr;
    _dpos = dinfo._entpos;
    _dpdir = dinfo._entpdir;
    _dtent = dinfo._tent;
    _dtenterr = dinfo._tenterr;
    _dentf = dinfo._entflt;
  }

  void
  Reflect::fillParentInfo(art::Ptr<SimParticle> sp) {
    GeomHandle<DetectorSystem> det;
    if(!sp.isNull()){
      _emom = sp->startMomentum().vect();
      _et0 = sp->startGlobalTime();
      _epos = det->toDetector(sp->startPosition());
// find the ultimate parent info
      art::Ptr<SimParticle> parsp = sp;
      while(!parsp->parent().isNull()){
	parsp = parsp->parent();
      }
      _ppdg = parsp->pdgId();
      CLHEP::Hep3Vector ppos = det->toDetector(parsp->startPosition());
      _ppos = ppos;
      _pt0 = sp->startGlobalTime();
      CLHEP::Hep3Vector pmom = parsp->startMomentum().vect();
      _pmom = pmom;
// project position into y=0 plane
      double flen = -ppos.y()/pmom.y();
      CLHEP::Hep3Vector planepos = ppos + flen*pmom;
      _pppos = planepos;
    }
  }

}  // end namespace md2e

// Part of the magic that makes this class a moddle.
// create an instance of the module.  It also registers
using mu2e::Reflect;
DEFINE_ART_MODULE(Reflect);
