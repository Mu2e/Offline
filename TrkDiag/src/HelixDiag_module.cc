//
// Helix Fit diagnostics
//
// Original author D. Brown and G. Tassielli
//

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
// mu2e
#include "GeneralUtilities/inc/Angles.hh"
#include "Mu2eUtilities/inc/MVATools.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
// diagnostics
#include "TrkDiag/inc/TrkMCTools.hh"
#include "TrkReco/inc/TrkUtilities.hh"
#include "RecoDataProducts/inc/HelixHit.hh"
#include "TrkDiag/inc/HelixHitInfo.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/HelixSeedCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "TrkReco/inc/TrkTimeCalculator.hh"
// root
#include "TGraph.h"
#include "TH2F.h"
#include "TF1.h"
#include "TList.h"
#include "TArc.h"
#include "TTree.h"
#include "TMarker.h"
// C++
#include <functional>
#include <algorithm>
#include <iostream>
using namespace std; 
using CLHEP::Hep3Vector;

namespace mu2e {

  class HelixDiag : public art::EDAnalyzer {
    public:
      explicit HelixDiag(fhicl::ParameterSet const& pset);
      virtual ~HelixDiag();
      virtual void beginRun(art::Run const& run) override;
      // This is called for each event.
      virtual void analyze(art::Event const& e) override;
    private:
// config parameters
      int _diag;
      StrawHitFlag _dontuseflag;
      bool _mcdiag;
      int _minnprimary; // minimum # Primary hits to make plots
      double _targetradius;
      bool _plot;
      TrkFitFlag _plotinc; // inclusive conditions for plotting
      TrkFitFlag _plotexc; // exclusive conditions for plotting
      double				_cradres; // average center resolution along center position (mm)
      double				_cperpres; // average center resolution perp to center position (mm)
      double				_radres; // average radial resolution for circle fit (mm)
      double				_phires; // average azimuthal resolution on circle (rad)
      // event object tags      art::InputTag _shTag;
      art::InputTag _shTag;
      art::InputTag _shpTag;
      art::InputTag _shfTag;
      art::InputTag _hsTag;
      art::InputTag _mcdigisTag;
      art::InputTag _vdmcstepsTag;
      // cache of event objects
      const StrawHitCollection* _shcol;
      const StrawHitPositionCollection* _shpcol;
      const StrawHitFlagCollection* _shfcol;
      const HelixSeedCollection* _hscol;
      const StrawDigiMCCollection* _mcdigis;
      const StepPointMCCollection* _vdmcsteps;
      // reco offsets
      TrkTimeCalculator			_ttcalc;
      // mc time offsets
      SimParticleTimeOffset _toff;
      // Virtual Detector IDs
      vector<int> _midvids;
      // cache of BField at 0,0,0
      double _bz0;
      // helper functions
      bool findData(const art::Event& e);
      void fillHitInfoMC(art::Ptr<SimParticle> const& pspp, StrawDigiMC const& digimc, HitInfoMC& hinfomc);
      bool fillMCHelix(art::Ptr<SimParticle> const& pspp);
      // display functions
      void plotXY(art::Ptr<SimParticle> const& pspp, HelixSeed const& myseed, unsigned ihel);
      void plotZ(art::Ptr<SimParticle> const& pspp, HelixSeed const& myseed, unsigned ihel);
      bool use(HelixHit const& hhit) const;
      bool stereo(HelixHit const& hhit) const;
      bool primary(art::Ptr<SimParticle> const& pspp, HelixHit const& hhit);
       // TTree and branch variables
      TTree *_hdiag;
      Int_t _iev;
      Bool_t _hitsOK, _circleInit, _phizInit, _circleOK, _phizOK, _helixOK, _mchelixOK;
      Float_t _mct0;
      Bool_t _circleConverged, _phizConverged, _helixConverged;
      Float_t _tct0, _tct0err, _ht0, _ht0err;
      RobustHelix _rhel;
      Int_t _nhits, _nused, _nprimary, _npused, _nptot;
      Int_t _pdg, _gen, _proc;
      vector<HelixHitInfo> _hhinfo;
      vector<HelixHitInfoMC> _hhinfomc;
      RobustHelix _mch;
      Float_t _mcmom, _mcpz;
  };

  HelixDiag::~HelixDiag() {
  }

  HelixDiag::HelixDiag(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _diag(pset.get<int>("DiagLevel",1)),
    _dontuseflag(pset.get<vector<string>>("UseFlag",vector<string>{"Outlier"})),
    _mcdiag(pset.get<bool>("MonteCarloDiag",true)),
    _minnprimary(pset.get<int>("MinimumPrimaryHits",10)),
    _targetradius(pset.get<double>("TargetRadius",75)),
    _plot(pset.get<bool>("PlotHelices",false)),
    _plotinc              (pset.get<vector<string> >("InclusivePlotFlagBits",vector<string>{"HitsOK"})),
    _plotexc              (pset.get<vector<string> >("ExclusivePlotFlagBits",vector<string>{})),
    _cradres	 (pset.get<double>("CenterRadialResolution",12.0)),
    _cperpres	 (pset.get<double>("CenterPerpResolution",12.0)),
    _radres	 (pset.get<double>("RadiusResolution",10.0)),
    _phires	 (pset.get<double>("AzimuthREsolution",0.1)),
    _shTag(pset.get<string>("StrawHitCollectionTag","makeSH")),
    _shpTag(pset.get<string>("StrawHitPositionCollectionTag","MakeStereoHits")),
    _shfTag(pset.get<string>("StrawHitFlagCollectionTag","PosHelixFinder")),
    _hsTag(pset.get<string>("HelixSeedCollectionTag","PosHelixFinder")),
    _mcdigisTag(pset.get<art::InputTag>("StrawDigiMCCollection","makeSH")),
    _vdmcstepsTag(pset.get<art::InputTag>("VDStepPointMCCollection","detectorFilter:virtualdetector")),
    _ttcalc            (pset.get<fhicl::ParameterSet>("T0Calculator",fhicl::ParameterSet())),
    _toff(pset.get<fhicl::ParameterSet>("TimeOffsets"))
  {
    if(_diag > 0){
      art::ServiceHandle<art::TFileService> tfs;
      _hdiag=tfs->make<TTree>("hdiag","Helix Finding diagnostics");
      _hdiag->Branch("iev",&_iev,"iev/I");
      _hdiag->Branch("hitsOK",&_hitsOK,"hitsOK/B");
      _hdiag->Branch("circleInit",&_circleInit,"circleInit/B");
      _hdiag->Branch("phizInit",&_phizInit,"phizInit/B");
      _hdiag->Branch("circleOK",&_circleOK,"circleOK/B");
      _hdiag->Branch("phizOK",&_phizOK,"phizOK/B");
      _hdiag->Branch("helixOK",&_helixOK,"helixOK/B");
      _hdiag->Branch("mchelixOK",&_mchelixOK,"mchelixOK/B");
      _hdiag->Branch("mct0",&_mct0,"mct0/F");
      _hdiag->Branch("circleConverged",&_circleConverged,"circleConverged/B");
      _hdiag->Branch("phizConverged",&_phizConverged,"phizConverged/B");
      _hdiag->Branch("helixConverged",&_helixConverged,"helixConverged/B");
      _hdiag->Branch("tct0",&_tct0,"tct0/F");
      _hdiag->Branch("tct0err",&_tct0err,"tct0err/F");
      _hdiag->Branch("ht0",&_ht0,"ht0/F");
      _hdiag->Branch("ht0err",&_ht0err,"ht0err/F");
      _hdiag->Branch("rhel",&_rhel);
      _hdiag->Branch("nhits",&_nhits,"nhits/I");
      _hdiag->Branch("nused",&_nused,"nused/I");
      if(_mcdiag){
	_hdiag->Branch("mch",&_mch);
	_hdiag->Branch("mcmom",&_mcmom,"mcmom/F");
	_hdiag->Branch("mcpz",&_mcpz,"mcpz/F");
	_hdiag->Branch("nprimary",&_nprimary,"nprimary/I");
	_hdiag->Branch("nptot",&_nptot,"nptot/I");
	_hdiag->Branch("npused",&_npused,"npused/I");
	_hdiag->Branch("pdg",&_pdg,"pdg/I");
	_hdiag->Branch("gen",&_gen,"gen/I");
	_hdiag->Branch("proc",&_proc,"proc/I");
      }
      if(_diag > 1){
	_hdiag->Branch("hh",&_hhinfo);
	if(_mcdiag){
	  _hdiag->Branch("hhmc",&_hhinfomc);
	}
      }
    }
  }

  void HelixDiag::beginRun(art::Run const& run){
  // multiple VIDs for the tracker midplane: can we ever fix this??
    _midvids.push_back(VirtualDetectorId::TT_Mid);
    _midvids.push_back(VirtualDetectorId::TT_MidInner);
    // get bfield
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    Hep3Vector vpoint_mu2e = det->toMu2e(Hep3Vector(0.0,0.0,0.0));
    _bz0 = bfmgr->getBField(vpoint_mu2e).z();
  }

  void HelixDiag::analyze(art::Event const& evt) {
    _iev=evt.id().event();
// find the data
    if(findData(evt)) {
     // loop over helices
      unsigned ihel(0);
      for(auto const& hseed : *_hscol) {
	// reset
	_hhinfo.clear();
	_hhinfomc.clear();
	_nprimary = _nptot = _npused = 0;

	RobustHelix const& rhel = hseed._helix;
	HelixHitCollection const& hhits = hseed._hhits;
	_nhits = hhits.size();
	TrkFitFlag const& status = hseed._status;
	_rhel = rhel;
	_tct0 = hseed.timeCluster()->t0().t0();
	_tct0err = hseed.timeCluster()->t0().t0Err();
	_ht0 = hseed._t0.t0();
	_ht0err = hseed._t0.t0Err();
	_hitsOK = status.hasAllProperties(TrkFitFlag::hitsOK);
	_circleInit = status.hasAllProperties(TrkFitFlag::circleInit);
	_phizInit = status.hasAllProperties(TrkFitFlag::phizInit);
	_circleOK = status.hasAllProperties(TrkFitFlag::circleOK);
	_phizOK = status.hasAllProperties(TrkFitFlag::phizOK);
	_helixOK = status.hasAllProperties(TrkFitFlag::helixOK);
	_circleConverged = status.hasAllProperties(TrkFitFlag::circleConverged);
	_phizConverged = status.hasAllProperties(TrkFitFlag::phizConverged);
	_helixConverged = status.hasAllProperties(TrkFitFlag::helixConverged);
	vector<StrawHitIndex> hits;
	art::Ptr<SimParticle> pspp;
	_nused = 0;
	for(auto const& hhit : hhits){
	  hits.push_back(hhit._shidx);
	  if(!hhit._flag.hasAnyProperty(StrawHitFlag::outlier))++_nused;
	}
	if(_mcdiag) {
	  // get information about the primary particle (produced most hits)
	  _nprimary = TrkMCTools::primaryParticle(pspp,hits,_mcdigis);
	  _nptot = TrkMCTools::countDigis(pspp,_mcdigis);
	  _pdg = pspp->pdgId();
	  _proc = pspp->realCreationCode();
	  if(pspp->genParticle().isNonnull())
	    _gen = pspp->genParticle()->generatorId().id();
	  else
	    _gen = -1;
	  // fill MC true helix parameters
	  _mchelixOK = fillMCHelix(pspp);
	  _mct0 = 0.0;
	  unsigned nmc = 0;
	  _npused = 0;
	  for(auto const& hhit : hhits) {
	    StrawDigiMC const& mcdigi = _mcdigis->at(hhit._shidx);
	    art::Ptr<StepPointMC> spmcp;
	    if (TrkMCTools::stepPoint(spmcp,mcdigi) >= 0 &&
	      spmcp->simParticle() == pspp ){
	      ++nmc;
	      _mct0 += _toff.timeWithOffsetsApplied(*spmcp);
	      if(!hhit._flag.hasAnyProperty(StrawHitFlag::outlier))++_npused;
	    }
	    if(_diag > 1){
	      HelixHitInfoMC hhinfomc;
	      fillHitInfoMC(pspp,mcdigi,hhinfomc);
	      Hep3Vector mchpos = hhit.pos(); // sets z position
	      _mch.position(mchpos);
	      hhinfomc._hpos = mchpos;
	      hhinfomc._hphi = _mch.circleAzimuth(hhit.pos().z());

	      Hep3Vector mcdh = hhit.pos() - mchpos;
	      hhinfomc._dwire = mcdh.dot(hhit.wdir());
	      static Hep3Vector zaxis(0.0,0.0,1.0); // unit in z direction
	      Hep3Vector wtdir = zaxis.cross(hhit.wdir()); // transverse direction to the wire
	      hhinfomc._dtrans = mcdh.dot(wtdir);

	      _hhinfomc.push_back(hhinfomc);
	    }
	  }
	  if(nmc > 0)_mct0 /= nmc;
	}
	// plot if requested and the fit satisfies the requirements
	if( _plot && hseed._status.hasAllProperties(_plotinc) &&
	    (!hseed._status.hasAnyProperty(_plotexc))  &&
	    _nprimary >= _minnprimary) {
	    // fill graphs for display
	    plotXY(pspp,hseed,ihel);
	    plotZ(pspp,hseed,ihel);
	}
	
	if(_diag > 1){
	  for(auto const& hhit : hhits) {
	    HelixHitInfo hhinfo;
	    hhinfo._outlier = hhit._flag.hasAnyProperty(StrawHitFlag::outlier);
	    hhinfo._stereo = hhit._flag.hasAnyProperty(StrawHitFlag::stereo);
	    hhinfo._tdiv = hhit._flag.hasAnyProperty(StrawHitFlag::tdiv);
	    hhinfo._delta = hhit._flag.hasAnyProperty(StrawHitFlag::delta);
	    hhinfo._esel = hhit._flag.hasAnyProperty(StrawHitFlag::energysel);
	    hhinfo._resphi = hhit._flag.hasAnyProperty(StrawHitFlag::resolvedphi);
	    hhinfo._hhphi = hhit._phi;
	    hhinfo._hhpos = hhit.pos();
	    hhinfo._werr = hhit.posRes(StrawHitPosition::wire);
	    hhinfo._terr = hhit.posRes(StrawHitPosition::trans);
	    hhinfo._dt = _shcol->at(hhit._shidx).time() - hseed._t0.t0() -_ttcalc.strawHitTimeOffset(hhit.pos().z());
	    Hep3Vector hpos = hhit.pos(); // this sets the z to the correct value
	    rhel.position(hpos);
	    hhinfo._hpos = hpos;
	    hhinfo._hphi = rhel.circleAzimuth(hhit.pos().z());
// compute the chisquared componentes for this hit
	    static Hep3Vector zaxis(0.0,0.0,1.0); // unit in z direction
	    Hep3Vector const& wdir = hhit.wdir();
	    Hep3Vector wtdir = zaxis.cross(wdir); // transverse direction to the wire
	    Hep3Vector cvec = (hhit.pos() - rhel.center()).perpPart(); // vector from the circle center to the hit
	    Hep3Vector cdir = cvec.unit(); // direction from the circle center to the hit
	    Hep3Vector cperp = zaxis.cross(cdir); // direction perp to the radius
	    hhinfo._whdot = wdir.dot(cperp); // compare wire and circle tangent directions
	    hhinfo._hrho = cvec.mag(); // radius of this hit WRT the circle center

	    // positions
	    Hep3Vector dh = hhit.pos() - hpos;
	    double dwire = dh.dot(wdir); // projection along wire direction
	    double dtrans = dh.dot(wtdir); // transverse projection

	    double wres2 = std::pow(hhit.posRes(StrawHitPosition::wire),(int)2) +
	      std::pow(_cradres*cdir.dot(wdir),(int)2) +
	      std::pow(_cperpres*cperp.dot(wdir),(int)2);
	    double wtres2 = std::pow(hhit.posRes(StrawHitPosition::trans),(int)2) +
	      std::pow(_cradres*cdir.dot(wtdir),(int)2) +
	      std::pow(_cperpres*cperp.dot(wtdir),(int)2);

	    hhinfo._dwire = dwire;
	    hhinfo._dtrans = dtrans;
	    hhinfo._wres = sqrt(wres2);
	    hhinfo._wtres = sqrt(wtres2);
	    hhinfo._chisq = sqrt( dwire*dwire/wres2 + dtrans*dtrans/wtres2 );
	    hhinfo._hqual = hhit._hqual;

	    _hhinfo.push_back(hhinfo);
	  }
	}
	  // fill the tree
	_hdiag->Fill();
	// if requested, plot the hits and helices
	++ihel;
      }
    } else
      cout << "HelixDiag_module can't find data" << endl;
  }

  bool HelixDiag::findData(const art::Event& evt){
    _shcol = 0;
    _shpcol = 0;
    _shfcol = 0;
    _hscol = 0;
    _mcdigis = 0;
    _vdmcsteps = 0;
// nb: getValidHandle does the protection (exception) on handle validity so I don't have to
    auto shH = evt.getValidHandle<StrawHitCollection>(_shTag);
    _shcol = shH.product();
    auto shpH = evt.getValidHandle<StrawHitPositionCollection>(_shpTag);
    _shpcol = shpH.product();
    auto shfH = evt.getValidHandle<StrawHitFlagCollection>(_shfTag);
    _shfcol = shfH.product();
    auto hsH = evt.getValidHandle<HelixSeedCollection>(_hsTag);
    _hscol = hsH.product();
    if(_mcdiag){
      auto mcdH = evt.getValidHandle<StrawDigiMCCollection>(_mcdigisTag);
      _mcdigis = mcdH.product();
      auto mcstepsH = evt.getValidHandle<StepPointMCCollection>(_vdmcstepsTag);
      _vdmcsteps = mcstepsH.product();
      // update time offsets
      _toff.updateMap(evt);
    }

    return _shcol != 0 && _shpcol != 0 && _shfcol != 0 && _hscol != 0 
      && ((_mcdigis != 0 && _vdmcsteps != 0 ) || !_mcdiag);
  }


  void HelixDiag::plotXY(art::Ptr<SimParticle> const& pspp, HelixSeed const& hseed, unsigned ihel) {
    RobustHelix const& rhel = hseed._helix;
    HelixHitCollection const& hhits = hseed._hhits;

    static unsigned igraph = 0;
    igraph++;
    art::ServiceHandle<art::TFileService> tfs;
    char ce_stereo_used_name[100];
    snprintf(ce_stereo_used_name,100,"ce_stereo_used_shxy%i",igraph);
    char ce_stereo_notused_name[100];
    snprintf(ce_stereo_notused_name,100,"ce_stereo_notused_shxy%i",igraph);
    char ce_notstereo_used_name[100];
    snprintf(ce_notstereo_used_name,100,"ce_notstereo_used_shxy%i",igraph);
    char ce_notstereo_notused_name[100];
    snprintf(ce_notstereo_notused_name,100,"ce_notstereo_notused_shxy%i",igraph);
    char bkg_stereo_used_name[100];
    snprintf(bkg_stereo_used_name,100,"bkg_stereo_used_shxy%i",igraph);
    char bkg_stereo_notused_name[100];
    snprintf(bkg_stereo_notused_name,100,"bkg_stereo_notused_shxy%i",igraph);
    char bkg_notstereo_used_name[100];
    snprintf(bkg_notstereo_used_name,100,"bkg_notstereo_used_shxy%i",igraph);
    char bkg_notstereo_notused_name[100];
    snprintf(bkg_notstereo_notused_name,100,"bkg_notstereo_notused_shxy%i",igraph);
    char title[100];
    snprintf(title,100,"StrawHit XY evt %i hel %i;mm;rad",_iev,ihel);
    TH2F* ce_stereo_used = tfs->make<TH2F>(ce_stereo_used_name,title,100,-500,500,100,-500,500);
    TH2F* ce_stereo_notused = tfs->make<TH2F>(ce_stereo_notused_name,title,100,-500,500,100,-500,500);
    TH2F* ce_notstereo_used = tfs->make<TH2F>(ce_notstereo_used_name,title,100,-500,500,100,-500,500);
    TH2F* ce_notstereo_notused = tfs->make<TH2F>(ce_notstereo_notused_name,title,100,-500,500,100,-500,500);
    TH2F* bkg_stereo_used = tfs->make<TH2F>(bkg_stereo_used_name,title,100,-500,500,100,-500,500);
    TH2F* bkg_stereo_notused = tfs->make<TH2F>(bkg_stereo_notused_name,title,100,-500,500,100,-500,500);
    TH2F* bkg_notstereo_used = tfs->make<TH2F>(bkg_notstereo_used_name,title,100,-500,500,100,-500,500);
    TH2F* bkg_notstereo_notused = tfs->make<TH2F>(bkg_notstereo_notused_name,title,100,-500,500,100,-500,500);

    ce_stereo_used->SetMarkerStyle(kFullTriangleUp);
    ce_stereo_used->SetMarkerColor(kRed);
    ce_stereo_notused->SetMarkerStyle(kOpenTriangleUp);
    ce_stereo_notused->SetMarkerColor(kRed);
    ce_notstereo_used->SetMarkerStyle(kFullCircle);
    ce_notstereo_used->SetMarkerColor(kRed);
    ce_notstereo_notused->SetMarkerStyle(kOpenCircle);
    ce_notstereo_notused->SetMarkerColor(kRed);
    bkg_stereo_used->SetMarkerStyle(kFullTriangleUp);
    bkg_stereo_used->SetMarkerColor(kGreen);
    bkg_stereo_notused->SetMarkerStyle(kOpenTriangleUp);
    bkg_stereo_notused->SetMarkerColor(kGreen);
    bkg_notstereo_used->SetMarkerStyle(kFullCircle);
    bkg_notstereo_used->SetMarkerColor(kGreen);
    bkg_notstereo_notused->SetMarkerStyle(kOpenCircle);
    bkg_notstereo_notused->SetMarkerColor(kGreen);

    for(auto const& hhit : hhits ) {
      if(primary(pspp,hhit)){
	if (use(hhit) ) {
	  if (stereo(hhit)) {
	    ce_stereo_used->Fill(hhit._pos.x()-rhel.centerx(),hhit._pos.y()-rhel.centery());
	  }
	  else {
	    ce_notstereo_used->Fill(hhit._pos.x()-rhel.centerx(),hhit._pos.y()-rhel.centery());
	  }	      
	}
	else {
	  if (stereo(hhit)) {
	    ce_stereo_notused->Fill(hhit._pos.x()-rhel.centerx(),hhit._pos.y()-rhel.centery());
	  }
	  else {
	    ce_notstereo_notused->Fill(hhit._pos.x()-rhel.centerx(),hhit._pos.y()-rhel.centery());
	  }
	}
      }
      else {
	if (use(hhit)) {
	  if (stereo(hhit)) {
	    bkg_stereo_used->Fill(hhit._pos.x()-rhel.centerx(),hhit._pos.y()-rhel.centery());
	  }
	  else {
	    bkg_notstereo_used->Fill(hhit._pos.x()-rhel.centerx(),hhit._pos.y()-rhel.centery());
	  }	      
	}
	else {
	  if (stereo(hhit)) {
	    bkg_stereo_notused->Fill(hhit._pos.x()-rhel.centerx(),hhit._pos.y()-rhel.centery());
	  }
	  else {
	    bkg_notstereo_notused->Fill(hhit._pos.x()-rhel.centerx(),hhit._pos.y()-rhel.centery());
	  }
	}
      }
    }

  // fit
    TArc* fitarc = new TArc(0.0,0.0,rhel.radius());
    fitarc->SetLineColor(kRed);
    fitarc->SetLineWidth(2);
    fitarc->SetFillStyle(0);
    // draw the detector boundaries
    static double innerrad(380.0);
    static double outerrad(680.0);
    TArc* indet = new TArc(-rhel.centerx(),-rhel.centery(),innerrad);
    TArc* outdet = new TArc(-rhel.centerx(),-rhel.centery(),outerrad);
    indet->SetLineColor(kBlue);
    indet->SetFillStyle(0);
    outdet->SetLineColor(kBlue);
    outdet->SetFillStyle(0);

    TArc* target = new TArc(-rhel.centerx(),-rhel.centery(),_targetradius);
    target->SetLineColor(kBlack);
    target->SetFillStyle(0);
    // add these to the plot
    TList* flist = ce_stereo_used->GetListOfFunctions();
    flist->Add(fitarc);
    flist->Add(indet);
    flist->Add(outdet);
    flist->Add(target);

    if (_mcdiag) {
      // mc truth
      TArc* mcarc = new TArc(_mch.centerx()-rhel.centerx(),_mch.centery()-rhel.centery(),_mch.radius());
      mcarc->SetLineColor(kBlue);
      mcarc->SetLineWidth(2);
      mcarc->SetFillStyle(0);
      flist->Add(mcarc);
      // Plot the MC true Primary hit positions
      char mctruth_name[100];
      snprintf(mctruth_name,100,"mctshxy%i",igraph);
      TH2F* mct = tfs->make<TH2F>(mctruth_name,title,100,-500,500,100,-500,500);
      mct->SetMarkerStyle(5);
      mct->SetMarkerColor(kMagenta);

      for(auto const& hhit : hhits ) {
	StrawDigiMC const& mcdigi = _mcdigis->at(hhit._shidx);
	art::Ptr<StepPointMC> spmcp;
	if (TrkMCTools::stepPoint(spmcp,mcdigi) >= 0 && spmcp->simParticle() == pspp)
	  mct->Fill(spmcp->position().x()-rhel.centerx(),spmcp->position().y()-rhel.centery());
      }
    }
  }

  void HelixDiag::plotZ(art::Ptr<SimParticle> const& pspp, HelixSeed const& hseed, unsigned ihel) {
    RobustHelix const& rhel = hseed._helix;
    HelixHitCollection const& hhits = hseed._hhits;

    static unsigned igraph = 0;
    igraph++;
    art::ServiceHandle<art::TFileService> tfs;

    char ce_stereo_used_name[100];
    snprintf(ce_stereo_used_name,100,"ce_stereo_used_shphiz%i",igraph);
    char ce_stereo_notused_name[100];
    snprintf(ce_stereo_notused_name,100,"ce_stereo_notused_shphiz%i",igraph);
    char ce_notstereo_used_name[100];
    snprintf(ce_notstereo_used_name,100,"ce_notstereo_used_shphiz%i",igraph);
    char ce_notstereo_notused_name[100];
    snprintf(ce_notstereo_notused_name,100,"ce_notstereo_notused_shphiz%i",igraph);
    char bkg_stereo_used_name[100];
    snprintf(bkg_stereo_used_name,100,"bkg_stereo_used_shphiz%i",igraph);
    char bkg_stereo_notused_name[100];
    snprintf(bkg_stereo_notused_name,100,"bkg_stereo_notused_shphiz%i",igraph);
    char bkg_notstereo_used_name[100];
    snprintf(bkg_notstereo_used_name,100,"bkg_notstereo_used_shphiz%i",igraph);
    char bkg_notstereo_notused_name[100];
    snprintf(bkg_notstereo_notused_name,100,"bkg_notstereo_notused_shphiz%i",igraph);
    char title[100];
    snprintf(title,100,"StrawHit #phi Z evt %i hel %i;mm;rad",_iev, ihel);
    TH2F* ce_stereo_used = tfs->make<TH2F>(ce_stereo_used_name,title,100,-1500,1500,100,-12.5,12.5);
    TH2F* ce_stereo_notused = tfs->make<TH2F>(ce_stereo_notused_name,title,100,-1500,1500,100,-12.5,12.5);
    TH2F* ce_notstereo_used = tfs->make<TH2F>(ce_notstereo_used_name,title,100,-1500,1500,100,-12.5,12.5);
    TH2F* ce_notstereo_notused = tfs->make<TH2F>(ce_notstereo_notused_name,title,100,-1500,1500,100,-12.5,12.5);
    TH2F* bkg_stereo_used = tfs->make<TH2F>(bkg_stereo_used_name,title,100,-1500,1500,100,-12.5,12.5);
    TH2F* bkg_stereo_notused = tfs->make<TH2F>(bkg_stereo_notused_name,title,100,-1500,1500,100,-12.5,12.5);
    TH2F* bkg_notstereo_used = tfs->make<TH2F>(bkg_notstereo_used_name,title,100,-1500,1500,100,-12.5,12.5);
    TH2F* bkg_notstereo_notused = tfs->make<TH2F>(bkg_notstereo_notused_name,title,100,-1500,1500,100,-12.5,12.5);

    ce_stereo_used->SetMarkerStyle(kFullTriangleUp);
    ce_stereo_used->SetMarkerColor(kRed);
    ce_stereo_notused->SetMarkerStyle(kOpenTriangleUp);
    ce_stereo_notused->SetMarkerColor(kRed);
    ce_notstereo_used->SetMarkerStyle(kFullCircle);
    ce_notstereo_used->SetMarkerColor(kRed);
    ce_notstereo_notused->SetMarkerStyle(kOpenCircle);
    ce_notstereo_notused->SetMarkerColor(kRed);
    bkg_stereo_used->SetMarkerStyle(kFullTriangleUp);
    bkg_stereo_used->SetMarkerColor(kGreen);
    bkg_stereo_notused->SetMarkerStyle(kOpenTriangleUp);
    bkg_stereo_notused->SetMarkerColor(kGreen);
    bkg_notstereo_used->SetMarkerStyle(kFullCircle);
    bkg_notstereo_used->SetMarkerColor(kGreen);
    bkg_notstereo_notused->SetMarkerStyle(kOpenCircle);
    bkg_notstereo_notused->SetMarkerColor(kGreen);

    for(auto const& hhit : hhits) {
      if(primary(pspp,hhit)){
	if (use(hhit)) {
	  if (stereo(hhit)) {
	    ce_stereo_used->Fill(hhit._pos.z(),hhit._phi);
	  }
	  else {
	    ce_notstereo_used->Fill(hhit._pos.z(),hhit._phi);
	  }	      
	}
	else {
	  if (stereo(hhit)) {
	    ce_stereo_notused->Fill(hhit._pos.z(),hhit._phi);
	  }
	  else {
	    ce_notstereo_notused->Fill(hhit._pos.z(),hhit._phi);
	  }
	}
      }
      else {
	if (use(hhit)) {
	  if (stereo(hhit)) {
	    bkg_stereo_used->Fill(hhit._pos.z(),hhit._phi);
	  }
	  else {
	    bkg_notstereo_used->Fill(hhit._pos.z(),hhit._phi);
	  }	      
	}
	else {
	  if (stereo(hhit)) {
	    bkg_stereo_notused->Fill(hhit._pos.z(),hhit._phi);
	  }
	  else {
	    bkg_notstereo_notused->Fill(hhit._pos.z(),hhit._phi);
	  }
	}
      }
    }
    TF1* fitline = new TF1("fitline","[0]+[1]*x",-1500,1500);
    fitline->SetParameter(0,rhel.fz0());
    fitline->SetParameter(1,1.0/rhel.lambda());
    fitline->SetLineColor(kRed);
    TList* flist = ce_stereo_used->GetListOfFunctions();
    flist->Add(fitline);

    if (_mcdiag) {
      TF1* mcline = new TF1("mcline","[0]+[1]*x",-1500,1500);
      mcline->SetParameter(0,_mch.fz0());
      mcline->SetParameter(1,1.0/_mch.lambda());
      mcline->SetLineColor(kBlue);
      TList* flist = ce_stereo_used->GetListOfFunctions();
      flist->Add(mcline);
      // Plot the MC true Primary hits
      char mctruth_name[100];
      snprintf(mctruth_name,100,"mctshphiz%i",igraph);
      TH2F* mct = tfs->make<TH2F>(mctruth_name,title,100,-1500,1500,100,-12.5,12.5);
      mct->SetMarkerStyle(5);
      mct->SetMarkerColor(kMagenta);

      for(auto& hhit : hhits ) {
	StrawDigiMC const& mcdigi = _mcdigis->at(hhit._shidx);
	art::Ptr<StepPointMC> spmcp;
	if (TrkMCTools::stepPoint(spmcp,mcdigi) >= 0 && spmcp->simParticle() == pspp)
	  mct->Fill(spmcp->position().z(),spmcp->position().phi());
      }
    }
  }

   void HelixDiag::fillHitInfoMC(const art::Ptr<SimParticle>& pspp, StrawDigiMC const& digimc, HitInfoMC& hinfomc) {
    hinfomc.reset();
    art::Ptr<SimParticle> spp;
    art::Ptr<StepPointMC> spmcp;
    if(TrkMCTools::simParticle(spp,digimc) > 0 && TrkMCTools::stepPoint(spmcp,digimc) >= 0 ){
      hinfomc._pdg = spp->pdgId();
      hinfomc._proc = spp->realCreationCode();
      if(spp->genParticle().isNonnull())
	hinfomc._gen = spp->genParticle()->generatorId().id();
      hinfomc._rel = MCRelationship::relationship(pspp,spp);
      hinfomc._t0 = _toff.timeWithOffsetsApplied(*spmcp);
    }
  }

  bool HelixDiag::fillMCHelix(art::Ptr<SimParticle> const& pspp) {
    bool retval(false);
    GeomHandle<DetectorSystem> det;
    GlobalConstantsHandle<ParticleDataTable> pdt;
  // find the earliest step associated with this particle passing the tracker midplane
    cet::map_vector_key trkid = pspp->id();
    auto jmc = _vdmcsteps->end();
    for(auto imc = _vdmcsteps->begin();imc != _vdmcsteps->end(); ++imc ) {
    // find matching steps
      if(  imc->trackId() == trkid &&
	  (find(_midvids.begin(),_midvids.end(),imc->volumeId()) != _midvids.end()) ) {
//	  cout << "Found matching step " << *imc << endl;
	  if(jmc == _vdmcsteps->end() || imc->time() < jmc->time())
	    jmc = imc;
      }
    }
    if(jmc != _vdmcsteps->end()){
      // get momentum and position from this point
      _mcmom = jmc->momentum().mag();
      _mcpz = jmc->momentum().z();
      Hep3Vector pos = det->toDetector(jmc->position());
      double charge = pdt->particle(pspp->pdgId()).ref().charge();
      TrkUtilities::RobustHelixFromMom(pos,jmc->momentum(),charge,_bz0,_mch);
      retval = true;
    }
    return retval;
  }

  bool
  HelixDiag::use(HelixHit const& hhit) const {
     return !hhit._flag.hasAnyProperty(_dontuseflag);
  }

  bool
  HelixDiag::stereo(HelixHit const& hhit) const {
    static StrawHitFlag stereo(StrawHitFlag::stereo);
    return hhit._flag.hasAllProperties(stereo);
  }

  bool
  HelixDiag::primary(art::Ptr<SimParticle> const& pspp, HelixHit const& hhit) {
    StrawDigiMC const& mcdigi = _mcdigis->at(hhit._shidx);
    art::Ptr<StepPointMC> spmcp;
    return TrkMCTools::stepPoint(spmcp,mcdigi) >= 0 && spmcp->simParticle() == pspp;
  }

}

using mu2e::HelixDiag;
DEFINE_ART_MODULE(HelixDiag);
