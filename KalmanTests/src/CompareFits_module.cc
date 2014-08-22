//_p// compare 2 particle fits of the same track
//
// $Id: CompareFits_module.cc,v 1.5 2014/08/22 19:55:50 brownd Exp $
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

  class CompareFits : public art::EDAnalyzer {
  public:

    explicit CompareFits(fhicl::ParameterSet const& pset);
    virtual ~CompareFits() { }

    void beginJob();
    void analyze(const art::Event& e);

  private:

    // Module label of the module that performed the fits.
    std::string _fitterPrefix;
    std::string _extModName;
    std::string _pdname, _pmname; // primary particle data product names
    std::string _sname, _smname; // secondary particle data product names
    TrkParticle::type _pprimary, _psecondary; // particle PDG code
    TrkFitDirection::FitDirection _fitdir;
    double _maxdmom, _maxdtd;
    double _maxdp0, _maxdt0;
    double _zent; // z position of the entrance of the tracker
    unsigned _eventid;
    bool _extrapolate;
    // diagnostic of Kalman fit
    KalDiag _kdiag;
    // TTree for studying  fits
    TTree* _fitcomp;
    // TTree branches
    Int_t _ppart, _spart, _pfitstat, _sfitstat;
    Int_t _pnactive, _snactive;
    Float_t _pchisq, _schisq, _pfitcon, _sfitcon;
    Float_t _pt0, _pt0err, _st0, _st0err;
    Float_t _pmom, _pmomerr, _smom, _smomerr;
    helixpar _ppar, _pperr, _spar, _sperr;
    threevec _ppos, _spos, _ppdir, _spdir;
    Float_t _ptent, _stent, _ptenterr, _stenterr;
    Float_t _pentf, _sentf;
    MCTrkInfo _pmcinfo,_smcinfo;
// create 
    void createTree();
// fill tree
    void fillTree(FitInfo const& pinfo, FitInfo const& sinfo);
// fill fit information
    void fillFitInfo(const KalRep* krep,FitInfo& fitinfo) const;
// Function to select 'the same' tracks
    bool sameFit(FitInfo const& pinfo, FitInfo const& sinfo) const;
    void getEntranceZ();
  };

  CompareFits::CompareFits(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset),
    _fitterPrefix(pset.get<string>("fitterPrefix","TPR")),
    _extModName(pset.get<string>("ExtrapolationModuleName","TrkExt")),
    _pprimary((TrkParticle::type)pset.get<int>("PrimaryParticleCode",TrkParticle::e_minus)),
    _psecondary((TrkParticle::type)pset.get<int>("SecondaryParticleCode",TrkParticle::mu_minus)),
    _fitdir((TrkFitDirection::FitDirection)pset.get<int>("FitDirection",TrkFitDirection::downstream)),
    _maxdmom(pset.get<double>("MaxDeltaMom",10.0)),
    _maxdtd(pset.get<double>("MaxDeltaTanDip",0.1)),
    _maxdp0(pset.get<double>("MaxDeltaPhi0",0.1)),
    _maxdt0(pset.get<double>("MaxDeltaT0",5)),
    _kdiag(pset.get<fhicl::ParameterSet>("KalDiag",fhicl::ParameterSet())),
    _fitcomp(0)
  {
// construct the data product instance names for particles 
    TrkFitDirection fdir(_fitdir);
    TrkParticle primary(_pprimary);
    TrkParticle secondary(_psecondary);
    _pdname = fdir.name() +primary.name();
    _pmname =  _fitterPrefix + _pdname;
    _sname = fdir.name() +secondary.name();
    _smname =  _fitterPrefix + _sname;
  }

  void CompareFits::beginJob( ){
// create the tree
    createTree();
// initialize counter
    _eventid = 0;
  }

  void CompareFits::getEntranceZ( ){
 // get the virtual detector at the tracker entrance and take its z position
    GeomHandle<VirtualDetector> vdg;
    GeomHandle<DetectorSystem> det;
    CLHEP::Hep3Vector entpos = det->toDetector(vdg->getGlobal(VirtualDetectorId::TT_FrontPA));
    _zent = entpos.z();
  } 

  // For each event, look at tracker hits and calorimeter hits.
  void CompareFits::analyze(const art::Event& event) {
    if(_eventid==0)getEntranceZ();
    _eventid++;
// get MC info
    bool hasmc = _kdiag.findMCData(event);
 
    art::Handle<KalRepCollection> ptrksHandle;
    event.getByLabel(_pmname,_pdname,ptrksHandle);
    if(!ptrksHandle.isValid())return;
    KalRepCollection const& ptrks = *ptrksHandle;
//
    art::Handle<KalRepCollection> strksHandle;
    event.getByLabel(_smname,_sname,strksHandle);
    if(!strksHandle.isValid())return;
    KalRepCollection const& strks = *strksHandle;
// loop over pairs, and see if any match.
    for ( size_t ip=0; ip < ptrks.size(); ++ip ){
      KalRep const* pkrep = ptrks.get(ip);
      if ( pkrep != 0 && pkrep->fitCurrent() ){
	FitInfo pinfo;
	fillFitInfo(pkrep,pinfo);
	for ( size_t id=0; id < strks.size(); ++id ){
	  KalRep const* skrep = strks.get(id);
	  if ( skrep != 0 && skrep->fitCurrent() ){
	    FitInfo sinfo;
	    fillFitInfo(skrep,sinfo);
	    if(sameFit(pinfo,sinfo)){
	      // fill the branches with this info
	      fillTree(pinfo,sinfo);
	      _pmcinfo = _smcinfo = MCTrkInfo();
	      // get MC info for the upstream and downstream tracks
	      if(hasmc){
		art::Ptr<SimParticle> pmcinfo, smcinfo;
		_kdiag.findMCTrk(pkrep,pmcinfo);
		_kdiag.findMCTrk(skrep,smcinfo);
		// use these to find the points where the true particle enters the tracker
		if(pmcinfo.isNonnull() && smcinfo.isNonnull() && 
		    pmcinfo == smcinfo){
		  std::vector<MCStepItr> steps;
		  _kdiag.findMCSteps(_kdiag.mcData()._mcvdsteps,pmcinfo->id(),_kdiag.VDids(KalDiag::trackerEnt),steps);
		  if(steps.size() == 2){
		    // These are sorted by time: first should be upstream, second down
		    _kdiag.fillMCTrkInfo(steps[0],_pmcinfo);
		    _kdiag.fillMCTrkInfo(steps[1],_smcinfo);
		  } else
		    std::cout << "Didn't find 2 steps" << std::endl;
		} else
		  std::cout << "MC info doesn't match " << std::endl;
	      } else
		std::cout << "No MC info " << std::endl;
	      _fitcomp->Fill();
	    }
	  }
	}
      }
    }
  }

// fill fit information from a KalRep
  void CompareFits::fillFitInfo(const KalRep* krep,FitInfo& fitinfo) const {
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
    CompareFits::sameFit(FitInfo const& pinfo, FitInfo const& sinfo) const {
      // check basic info
      bool retval = pinfo._fitstatus > 0 && sinfo._fitstatus > 0;
      if(retval){
	// rough cuts on momentum and direction at tracker entrance
	double dmom = pinfo._entmom - sinfo._entmom;
	double dtd = pinfo._entpar._td - sinfo._entpar._td;
	double dp0 = pinfo._entpar._p0 - sinfo._entpar._p0;
	if(fabs(dp0) > M_PI){
	  if(dp0 > 0.0)
	    dp0 -= 2*M_PI;
	  else
	    dp0 += 2*M_PI;
	}
	double dt = pinfo._t0 - sinfo._t0;
	retval &= fabs(dmom) < _maxdmom && fabs(dtd) < _maxdtd &&
	  fabs(dp0) < _maxdp0 && fabs(dt) < _maxdt0;
      }
      return retval;
    }

  void
    CompareFits::createTree() {
      art::ServiceHandle<art::TFileService> tfs;
      _fitcomp=tfs->make<TTree>("compare","Compare Kalman fits");
      // create the branches
      _fitcomp->Branch("eventid",&_eventid,"eventid/I");
      // primary information
      _fitcomp->Branch("ppart",&_ppart,"ppart/I");
      _fitcomp->Branch("pfitstat",&_pfitstat,"pfitstat/I");
      _fitcomp->Branch("pnactive",&_pnactive,"pnactive/I");
      _fitcomp->Branch("pchisq",&_pchisq,"pchisq/F");
      _fitcomp->Branch("pfitcon",&_pfitcon,"pfitcon/F");
      _fitcomp->Branch("pt0",&_pt0,"pt0/F");
      _fitcomp->Branch("pt0err",&_pt0err,"pt0err/F");
      _fitcomp->Branch("pmom",&_pmom,"pmom/F");
      _fitcomp->Branch("pmomerr",&_pmomerr,"pmomerr/F");
      _fitcomp->Branch("ppar",&_ppar,"pd0/F:pp0/F:pom/F:pz0/F:ptd/F");
      _fitcomp->Branch("pperr",&_pperr,"pd0err/F:pp0err/F:pomerr/F:pz0err/F:ptderr/F");
      _fitcomp->Branch("ppos",&_ppos,"px/F:py/F:pz/F");
      _fitcomp->Branch("ppdir",&_ppdir,"ppx/F:ppy/F:ppz/F");
      _fitcomp->Branch("ptent",&_ptent,"ptent/F");
      _fitcomp->Branch("ptenterr",&_ptenterr,"ptenterr/F");
      _fitcomp->Branch("pentf",&_pentf,"pentf/F");
      _fitcomp->Branch("pmcinfo",&_pmcinfo,"pmcpdgid/I:pmctime/F:pmcmom/F:pmcx/F:pmcy/F:pmcz/F:pmcd0/F:pmcp0/F:pmcom/F:pmcz0/F:pmctd/F");
      // secondary information
      _fitcomp->Branch("spart",&_spart,"spart/I");
      _fitcomp->Branch("sfitstat",&_sfitstat,"sfitstat/I");
      _fitcomp->Branch("snactive",&_snactive,"snactive/I");
      _fitcomp->Branch("schisq",&_schisq,"schisq/F");
      _fitcomp->Branch("sfitcon",&_sfitcon,"sfitcon/F");
      _fitcomp->Branch("st0",&_st0,"st0/F");
      _fitcomp->Branch("st0err",&_st0err,"st0err/F");
      _fitcomp->Branch("smom",&_smom,"smom/F");
      _fitcomp->Branch("smomerr",&_smomerr,"smomerr/F");
      _fitcomp->Branch("spar",&_spar,"sd0/F:sp0/F:som/F:sz0/F:std/F");
      _fitcomp->Branch("sperr",&_sperr,"sd0err/F:sp0err/F:somerr/F:sz0err/F:stderr/F");
      _fitcomp->Branch("spos",&_spos,"sx/F:sy/F:sz/F");
      _fitcomp->Branch("spdir",&_spdir,"spx/F:spy/F:spz/F");
      _fitcomp->Branch("stent",&_stent,"stent/F");
      _fitcomp->Branch("stent",&_stent,"stent/F");
      _fitcomp->Branch("stenterr",&_stenterr,"stenterr/F");
      _fitcomp->Branch("sentf",&_sentf,"sentf/F");
      _fitcomp->Branch("smcinfo",&_smcinfo,"smcpdgid/I:smctime/F:smcmom/F:smcx/F:smcy/F:smcz/F:smcd0/F:smcp0/F:smcom/F:smcz0/F:smctd/F");
    }

  void
    CompareFits::fillTree(FitInfo const& pinfo, FitInfo const& sinfo) {
      _ppart = pinfo._fitpart;
      _pfitstat = pinfo._fitstatus;
      _pnactive = pinfo._nactive;
      _pchisq = pinfo._chisq;
      _pfitcon = pinfo._fitcon;
      _pt0 = pinfo._t0;
      _pt0err = pinfo._t0err;
      _pmom = pinfo._entmom;
      _pmomerr = pinfo._entmomerr;
      _ppar = pinfo._entpar;
      _pperr = pinfo._entperr;
      _ppos = pinfo._entpos;
      _ppdir = pinfo._entpdir;
      _ptent = pinfo._tent;
      _ptenterr = pinfo._tenterr;
      _pentf = pinfo._entflt;
      // downstream info
      _spart = sinfo._fitpart;
      _sfitstat = sinfo._fitstatus;
      _snactive = sinfo._nactive;
      _schisq = sinfo._chisq;
      _sfitcon = sinfo._fitcon;
      _st0 = sinfo._t0;
      _st0err = sinfo._t0err;
      _smom = sinfo._entmom;
      _smomerr = sinfo._entmomerr;
      _spar = sinfo._entpar;
      _sperr = sinfo._entperr;
      _spos = sinfo._entpos;
      _spdir = sinfo._entpdir;
      _stent = sinfo._tent;
      _stenterr = sinfo._tenterr;
      _sentf = sinfo._entflt;
    }

}  // end namespace md2e

// Part of the magic that makes this class a moddle.
// create an instance of the module.  It also registers
using mu2e::CompareFits;
DEFINE_ART_MODULE(CompareFits);
