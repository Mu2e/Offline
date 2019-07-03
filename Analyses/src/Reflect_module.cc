//
// Look for particles coming from the calorimeter and reflecting back in the
// magnetic mirror
//
// $Id: Reflect_module.cc,v 1.12 2014/09/10 18:49:17 brownd Exp $
// $Author: brownd $
// $Date: 2014/09/10 18:49:17 $
//
// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
// services
#include "BFieldGeom/inc/BFieldConfig.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
#include "GeometryService/inc/DetectorSystem.hh"
// data
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "RecoDataProducts/inc/TrkExtTrajCollection.hh"
#include "DataProducts/inc/XYZVec.hh"
// ROOT incldues
#include "TTree.h"
#include "TH1F.h"
// Need this for the BaBar headers.
using namespace CLHEP;
// BaBar includes
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
#include "BTrk/TrkBase/TrkHelixUtils.hh"
// mu2e tracking
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "TrkDiag/inc/KalDiag.hh"
#include "TrkDiag/inc/TrkInfo.hh"
// C++ includes.
#include <iostream>
#include <string>
#include <math.h>
#include <memory>
// This is fragile and needs to be last until CLHEP is
// properly qualified and included in the BaBar classes.
#include "RecoDataProducts/inc/KalRepCollection.hh"

using namespace std;

namespace mu2e {
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
    Int_t _trkid,_eventid, _runid, _subrunid;
    bool _extrapolate;
    // diagnostic of Kalman fits
    KalDiag _kdiag;
    // TTree for studying reflecting fits
    TTree* _reflect;
    // TTree branches
    TrkInfo _utrkinfo, _dtrkinfo;
    TrkFitInfo _utrkinfo_ent, _dtrkinfo_ent;
    TrkInfoMC _mcinfo;
    TrkInfoMCStep _umcinfo,_dmcinfo;
    XYZVec _opos, _omom, _ppos, _pmom, _pppos;
    Int_t _ppdg;
    Float_t _ot0, _pt0;
    Int_t _uextnpa,_dextnpa,_uextnst,_dextnst;
    Float_t _uextdppa,_dextdppa,_uextdpst,_dextdpst;
// create 
    void createTree();
    // Function to pair upstream and downstream fits
    bool reflection() const;
    void fillMCInfo(art::Ptr<SimParticle> sp);
  };

  Reflect::Reflect(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset),
    _fitterPrefix(pset.get<string>("fitterPrefix","TPR")),
    _extModName(pset.get<string>("ExtrapolationModuleName","TrkExt")),
    _ipart(pset.get<int>("ParticleCode",11)),
    _mindmom(pset.get<double>("MinDeltaMom",-20.0)),
    _maxdmom(pset.get<double>("MaxDeltaMom",20.0)),
    _maxdtd(pset.get<double>("MaxDeltaTanDip",0.5)),
    _maxdp0(pset.get<double>("MaxDeltaPhi0",0.5)),
    _mindt0(pset.get<double>("MinDeltaT0",-50)),
    _maxdt0(pset.get<double>("MaxDeltaT0",-500)),
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
  }

  // For each event, look at tracker hits and calorimeter hits.
  void Reflect::analyze(const art::Event& event) {
    _eventid = event.event();
    _runid = event.run();
    _subrunid = event.subRun();
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
	  _kdiag.fillTrkInfo(ukrep,_utrkinfo);
	  _kdiag.fillTrkFitInfo(ukrep,_utrkinfo_ent);
	  for ( size_t ide=0; ide < dtrks.size(); ++ide ){
	    KalRep const* dkrep = dtrks.get(ide);
	    if ( dkrep != 0 && dkrep->fitCurrent() ){
// fill fit info
	      _kdiag.fillTrkInfo(dkrep,_dtrkinfo);
	      _kdiag.fillTrkFitInfo(dkrep,_dtrkinfo_ent);
	      if(reflection()){
// fill the branches with this info
		_umcinfo.reset(); _dmcinfo.reset();
// get MC info for the upstream and downstream tracks
		if(hasmc){
		  art::Ptr<SimParticle> umcsp, dmcsp;
		  _kdiag.findMCTrk(ukrep,umcsp);
		  _kdiag.findMCTrk(dkrep,dmcsp);
// use these to find the points where the true particle enters the tracker
		  if(umcsp.isNonnull() && dmcsp.isNonnull() && 
		    umcsp == dmcsp){
		    // fill generic MC track information
		    _kdiag.fillTrkInfoMC(umcsp,0,_mcinfo);
		    std::vector<MCStepItr> steps;
		    _kdiag.findMCSteps(_kdiag.mcData()._mcvdsteps,umcsp->id(),_kdiag.VDids(KalDiag::trackerEnt),steps);
		    // fill specific MC information
		    fillMCInfo(umcsp);
		    if(steps.size() == 2){
// These are sorted by time: first should be upstream, second down
		      _kdiag.fillTrkInfoMCStep(steps[0],_umcinfo);
		      _kdiag.fillTrkInfoMCStep(steps[1],_dmcinfo);
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

  bool 
  Reflect::reflection() const {
// check basic info
    bool retval = _utrkinfo._status > 0 && _dtrkinfo._status > 0 &&
      _utrkinfo._pdg == _dtrkinfo._pdg;
    if(retval){
      // rough cuts on momentum and direction at tracker entrance
      double dmom = _utrkinfo_ent._fitmom-_dtrkinfo_ent._fitmom;
      double dtd = _utrkinfo_ent._fitpar._td + _dtrkinfo_ent._fitpar._td; // note sign change!!!
      double dp0 = _utrkinfo_ent._fitpar._p0 - _dtrkinfo_ent._fitpar._p0;
      if(fabs(dp0) > M_PI){
	if(dp0 > 0.0)
	  dp0 -= 2*M_PI;
	else
	  dp0 += 2*M_PI;
      }
      double dt = _utrkinfo._t0 - _dtrkinfo._t0;
      retval &= dmom > _mindmom && dmom < _maxdmom && fabs(dtd) < _maxdtd &&
	fabs(dp0) < _maxdp0 && dt < _mindt0 && dt > _maxdt0;
    }
    return retval;
  }


  void
  Reflect::fillMCInfo(art::Ptr<SimParticle> sp) {
    GeomHandle<DetectorSystem> det;
    if(!sp.isNull()){
      CLHEP::Hep3Vector opos = det->toDetector(sp->startPosition());
      CLHEP::Hep3Vector omom = sp->startMomentum().vect();
      _opos = opos;
      _omom = omom;
      _ot0 = sp->startGlobalTime();
// find the ultimate parent info
      art::Ptr<SimParticle> parsp = sp;
      while(!parsp->parent().isNull()){
       parsp = parsp->parent();
      }
      _ppdg = parsp->pdgId();
      CLHEP::Hep3Vector ppos = det->toDetector(parsp->startPosition());
      CLHEP::Hep3Vector pmom = parsp->startMomentum().vect();
      _ppos = ppos;
      _pmom = pmom;
      _pt0 = parsp->startGlobalTime();
// project position into y=0 plane
      double flen = -ppos.y()/pmom.y();
      CLHEP::Hep3Vector planepos = ppos + flen*pmom;
      _pppos = planepos;
    }
  }

  void
  Reflect::createTree() {
    art::ServiceHandle<art::TFileService> tfs;
    _reflect=tfs->make<TTree>("reflect","reflection diagnostics");
// create the branches
    _reflect->Branch("eventid",&_eventid,"eventid/I");
    _reflect->Branch("runid",&_runid,"runid/I");
    _reflect->Branch("subrunid",&_subrunid,"subrunid/I");
    // upstream and downstream track information
    _reflect->Branch("utrk",&_utrkinfo,TrkInfo::leafnames().c_str());
    _reflect->Branch("utrkent",&_utrkinfo_ent,TrkFitInfo::leafnames().c_str());
    _reflect->Branch("dtrk",&_dtrkinfo,TrkInfo::leafnames().c_str());
    _reflect->Branch("dtrkent",&_dtrkinfo_ent,TrkFitInfo::leafnames().c_str());
    _reflect->Branch("mc",&_mcinfo,TrkInfoMC::leafnames().c_str());
    _reflect->Branch("umc",&_umcinfo,TrkInfoMCStep::leafnames().c_str());
    _reflect->Branch("dmc",&_dmcinfo,TrkInfoMCStep::leafnames().c_str());
    _reflect->Branch("opos",&_opos,Geom::XYZnames().c_str());
    _reflect->Branch("omom",&_omom,Geom::XYZnames().c_str());
    _reflect->Branch("ot0",&_ot0,"ot0/F");
    _reflect->Branch("ppdg",&_ppdg,"ppdg/I");
    _reflect->Branch("ppos",&_ppos,Geom::XYZnames().c_str());
    _reflect->Branch("pmom",&_pmom,Geom::XYZnames().c_str());
    _reflect->Branch("pt0",&_pt0,"pt0/F");
    _reflect->Branch("pppos",&_pppos,Geom::XYZnames().c_str());
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
}  // end namespace md2e

// Part of the magic that makes this class a moddle.
// create an instance of the module.  It also registers
using mu2e::Reflect;
DEFINE_ART_MODULE(Reflect);
