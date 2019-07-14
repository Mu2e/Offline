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
#include "art_root_io/TFileService.h"
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
#include "TrkDiag/inc/HelixHitInfo.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
// data
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
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
#include "TEllipse.h"
#include "TMath.h"
#include "TMarker.h"
#include "Math/VectorUtil.h"
// C++
#include <functional>
#include <algorithm>
#include <iostream>
using namespace std; 
using CLHEP::Hep3Vector;
using namespace ROOT::Math::VectorUtil;

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
      bool _mcdiag, _mcsel, _useshfcol;
      StrawHitFlag  _hsel, _hbkg;
      unsigned _minnprimary; // minimum # Primary hits to make plots
      int _mcgen; // MC generator code of primary
      double _targetradius;
      bool _plot;
      bool _hcenter;
      double _xysize, _zsize, _fsize;
      TrkFitFlag _plotinc; // inclusive conditions for plotting
      TrkFitFlag _plotexc; // exclusive conditions for plotting
      double				_cradres; // average center resolution along center position (mm)
      double				_cperpres; // average center resolution perp to center position (mm)
      // event object tags      art::InputTag _shTag;
      art::InputTag _chTag;
      art::InputTag _hsTag;
      art::InputTag   _shfTag;
      art::InputTag _mcdigisTag;
      art::InputTag _vdmcstepsTag;
      art::InputTag _primaryTag;
      // cache of event objects
      const ComboHitCollection* _chcol;
      const HelixSeedCollection* _hscol;
      const StrawHitFlagCollection*	  _shfcol;
      const StrawDigiMCCollection* _mcdigis;
      const StepPointMCCollection* _vdmcsteps;
      const PrimaryParticle* _primary;
      // reco offsets
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
      void plotXY(const art::Event& e, art::Ptr<SimParticle> const& pspp, HelixSeed const& myseed, unsigned ihel);
      void plotZ(const art::Event& e, art::Ptr<SimParticle> const& pspp, HelixSeed const& myseed, unsigned ihel);
      bool use(ComboHit const& hhit) const;
      bool stereo(ComboHit const& hhit) const;
      bool primary(art::Ptr<SimParticle> const& pspp,size_t index);
      bool selectedHit(size_t index);
       // TTree and branch variables
      TTree *_hdiag;
      Int_t _iev;
      Bool_t _hitsOK, _circleInit, _phizInit, _circleOK, _phizOK, _helixOK, _mchelixOK;
      Float_t _mct0;
      Bool_t _circleConverged, _phizConverged, _helixConverged;
      Float_t _tct0, _tct0err, _ht0, _ht0err;
      RobustHelix _rhel;
      Int_t _nhits, _nused, _npused, _nptot;
      unsigned _nprimary;
      Int_t _pdg, _gen, _proc, _prel;
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
    _dontuseflag(pset.get<vector<string>>("DontUseFlag",vector<string>{"Outlier"})),
    _mcdiag(pset.get<bool>("MonteCarloDiag",true)),
    _mcsel(pset.get<bool>("MonteCarloSelection",true)),
    _useshfcol		(pset.get<bool>("UseFlagCollection")),
    _hsel(pset.get<std::vector<std::string> >("HitSelectionBits",vector<string>{"EnergySelection","TimeSelection","RadiusSelection"})),
    _hbkg(pset.get<vector<string> >("HitBackgroundBits",vector<string>{"Background"})),
    _minnprimary(pset.get<unsigned>("MinimumPrimaryHits",8)),
    _mcgen(pset.get<int>("MCGeneratorCode",2)),
    _targetradius(pset.get<double>("TargetRadius",75)),
    _plot(pset.get<bool>("PlotHelices",false)),
    _hcenter(pset.get<bool>("CenterPlotOnHelix",true)),
    _xysize(pset.get<double>("PlotXYSize",300.0)),
    _zsize(pset.get<double>("PlotZSize",1500.0)),
    _fsize(pset.get<double>("PlotPhiSize",2.0)),
    _plotinc              (pset.get<vector<string> >("InclusivePlotFlagBits",vector<string>{})),
    _plotexc              (pset.get<vector<string> >("ExclusivePlotFlagBits",vector<string>{})),
    _cradres	 (pset.get<double>("CenterRadialResolution",20.0)),
    _cperpres	 (pset.get<double>("CenterPerpResolution",12.0)),
    _chTag(pset.get<string>("ComboHitCollection")),
    _hsTag(pset.get<string>("HelixSeedCollection","HelixFinder:Positive")),
    _shfTag		(pset.get<string>("StrawHitFlagCollection")),
    _mcdigisTag(pset.get<art::InputTag>("StrawDigiMCCollection","makeSD")),
    _vdmcstepsTag(pset.get<art::InputTag>("VDStepPointMCCollection","detectorFilter:virtualdetector")),
    _primaryTag(pset.get<art::InputTag>("PrimaryParticleTag","FindMCPrimary")),
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
	_hdiag->Branch("mchelixOK",&_mchelixOK,"mchelixOK/B");
	_hdiag->Branch("mct0",&_mct0,"mct0/F");
	_hdiag->Branch("mch",&_mch);
	_hdiag->Branch("mcmom",&_mcmom,"mcmom/F");
	_hdiag->Branch("mcpz",&_mcpz,"mcpz/F");
	_hdiag->Branch("nprimary",&_nprimary,"nprimary/i");
	_hdiag->Branch("nptot",&_nptot,"nptot/I");
	_hdiag->Branch("npused",&_npused,"npused/I");
	_hdiag->Branch("pdg",&_pdg,"pdg/I");
	_hdiag->Branch("gen",&_gen,"gen/I");
	_hdiag->Branch("proc",&_proc,"proc/I");
	_hdiag->Branch("prel",&_prel,"prel/I");
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
	ComboHitCollection const& hhits = hseed._hhits;
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
	_nused = 0;
	std::vector<StrawDigiIndex> sdis;
	for(size_t ihh = 0;ihh < hhits.size(); ++ihh) {
	  ComboHit const& hhit = hhits[ihh];
	  hhits.fillStrawDigiIndices(evt,ihh,sdis);
	  if(!hhit.flag().hasAnyProperty(StrawHitFlag::outlier))_nused += hhit.nStrawHits();
	}
	art::Ptr<SimParticle> pspp;
	if(_mcdiag) {
	  _nprimary = 0; _pdg = -1; _proc = -1; _gen = -1; _prel = -1;
	  // get information about the primary particle (produced most hits)
	  MCRelationship mcrel;
	  TrkMCTools::primaryRelation(*_primary, *_mcdigis,sdis,
	  pspp, _nprimary, mcrel);
	  if(_nprimary >= _minnprimary){
	    _pdg = pspp->pdgId();
	    _proc = pspp->originParticle().creationCode();
	    _gen = _primary->primary().generatorId().id();
	    _nptot = TrkMCTools::countDigis(pspp,_mcdigis);
	    _prel = mcrel.relationship();
	    // fill MC true helix parameters
	    _mchelixOK = fillMCHelix(pspp);
	    _mct0 = 0.0;
	    unsigned nmc = 0;
	    _npused = 0;
	    for(size_t ihh = 0;ihh < hhits.size(); ++ihh) {
	      ComboHit const& hhit = hhits[ihh];
	      vector<StrawDigiIndex> sdis;
	      hhits.fillStrawDigiIndices(evt,ihh,sdis);
	      for(auto idigi : sdis) {
		StrawDigiMC const& mcdigi = _mcdigis->at(idigi);
		art::Ptr<StepPointMC> spmcp;
		if (TrkMCTools::stepPoint(spmcp,mcdigi) >= 0 &&
		    spmcp->simParticle() == pspp ){
		  ++nmc;
		  _mct0 += _toff.timeWithOffsetsApplied(*spmcp);
		  if(!hhit._flag.hasAnyProperty(StrawHitFlag::outlier))++_npused;
		}
	      }
	      if(_diag > 1){
		HelixHitInfoMC hhinfomc;
		StrawDigiMC const& mcdigi = _mcdigis->at(sdis[0]);
		fillHitInfoMC(pspp,mcdigi,hhinfomc);
		XYZVec mchpos = hhit.pos(); // sets z position
		_mch.position(mchpos);
		hhinfomc._hpos = mchpos;
		hhinfomc._hphi = _mch.circleAzimuth(hhit.pos().z());

		XYZVec mcdh = hhit.pos() - mchpos;
		hhinfomc._dwire = mcdh.Dot(hhit.wdir());
		static XYZVec zaxis(0.0,0.0,1.0); // unit in z direction
		XYZVec wtdir = zaxis.Cross(hhit.wdir()); // transverse direction to the wire
		hhinfomc._dtrans = mcdh.Dot(wtdir);
		_hhinfomc.push_back(hhinfomc);
	      }
	    }
	    if(nmc > 0)_mct0 /= nmc;
	  }
	}
	// plot if requested and the fit satisfies the requirements
	if( _plot && hseed._status.hasAllProperties(_plotinc) &&
	    (!hseed._status.hasAnyProperty(_plotexc))  &&
	    ((!_mcsel) || (_nprimary >= _minnprimary && _mcgen == _gen)) ) {
	  // fill graphs for display
	  plotXY(evt,pspp,hseed,ihel);
	  plotZ(evt,pspp,hseed,ihel);
	}

	if(_diag > 1){
	  for(auto const& hhit : hhits) {
	    HelixHitInfo hhinfo;
	    hhinfo._outlier = hhit.flag().hasAnyProperty(StrawHitFlag::outlier);
	    hhinfo._stereo = hhit.flag().hasAnyProperty(StrawHitFlag::stereo);
	    hhinfo._tdiv = hhit.flag().hasAnyProperty(StrawHitFlag::tdiv);
	    hhinfo._delta = hhit.flag().hasAnyProperty(StrawHitFlag::bkg);
	    hhinfo._esel = hhit.flag().hasAnyProperty(StrawHitFlag::energysel);
	    hhinfo._resphi = hhit.flag().hasAnyProperty(StrawHitFlag::resolvedphi);
	    hhinfo._hhphi = hhit.helixPhi();
	    hhinfo._hhpos = hhit.pos();
	    hhinfo._werr = hhit.posRes(ComboHit::wire);
	    hhinfo._terr = hhit.posRes(ComboHit::trans);
	    XYZVec hpos = hhit.pos(); // this sets the z to the correct value
	    rhel.position(hpos);
	    hhinfo._hpos = hpos;
	    hhinfo._hphi = rhel.circleAzimuth(hhit.pos().z());
	    // compute the chisquared componentes for this hit
	    XYZVec const& wdir = hhit.wdir();
	    XYZVec wtdir = Geom::ZDir().Cross(wdir); // transverse direction to the wire
	    XYZVec cvec = PerpVector(hhit.pos() - rhel.center(),Geom::ZDir()); // vector from the circle center to the hit
	    XYZVec cdir = cvec.unit(); // direction from the circle center to the hit
	    XYZVec cperp = Geom::ZDir().Cross(cdir); // direction perp to the radius
	    hhinfo._whdot = wdir.Dot(cdir); // compare wire and circle radius direction
	    hhinfo._hrho = sqrt(cvec.Mag2()); // radius of this hit WRT the circle center

	    // positions
	    XYZVec dh = hhit.pos() - hpos;
	    double dwire = dh.Dot(wdir); // projection along wire direction
	    double dtrans = dh.Dot(wtdir); // transverse projection

	    double wres2 = std::pow(hhit.posRes(ComboHit::wire),(int)2) +
	      std::pow(_cradres*cdir.Dot(wdir),(int)2) +
	      std::pow(_cperpres*cperp.Dot(wdir),(int)2);
	    double wtres2 = std::pow(hhit.posRes(ComboHit::trans),(int)2) +
	      std::pow(_cradres*cdir.Dot(wtdir),(int)2) +
	      std::pow(_cperpres*cperp.Dot(wtdir),(int)2);

	    hhinfo._dwire = dwire;
	    hhinfo._dtrans = dtrans;
	    hhinfo._wres = sqrt(wres2);
	    hhinfo._wtres = sqrt(wtres2);
	    hhinfo._chisq = dwire*dwire/wres2 + dtrans*dtrans/wtres2;
	    hhinfo._hqual = hhit.qual();

	    _hhinfo.push_back(hhinfo);
	  }
	}
	  // fill the tree
	if( (!_mcsel) || ( _nprimary >= _minnprimary && _mcgen == _gen))
	  _hdiag->Fill();
	++ihel;
      }
    } else
      cout << "HelixDiag_module can't find data" << endl;
  }

  bool HelixDiag::findData(const art::Event& evt){
    _chcol = 0;
    _shfcol = 0;
    _hscol = 0;
    _mcdigis = 0;
    _vdmcsteps = 0;
    _primary = 0;
// nb: getValidHandle does the protection (exception) on handle validity so I don't have to
    auto chH = evt.getValidHandle<ComboHitCollection>(_chTag);
    _chcol = chH.product();
    auto hsH = evt.getValidHandle<HelixSeedCollection>(_hsTag);
    _hscol = hsH.product();
    if(_mcdiag){
      auto mcdH = evt.getValidHandle<StrawDigiMCCollection>(_mcdigisTag);
      _mcdigis = mcdH.product();
      auto mcstepsH = evt.getValidHandle<StepPointMCCollection>(_vdmcstepsTag);
      _vdmcsteps = mcstepsH.product();
      // update time offsets
      _toff.updateMap(evt);
      auto pph = evt.getValidHandle<PrimaryParticle>(_primaryTag);
      _primary = pph.product();
    }
    if(_useshfcol){
      auto shfH = evt.getValidHandle<StrawHitFlagCollection>(_shfTag);
      _shfcol = shfH.product();
    }

    return _chcol != 0 && _hscol != 0 
      && ((_mcdigis != 0 && _vdmcsteps != 0 ) || !_mcdiag)
      && (_shfcol != 0 || (!_useshfcol));
  }


  void HelixDiag::plotZ(art::Event const& evt, art::Ptr<SimParticle> const& pspp, HelixSeed const& hseed, unsigned ihel) {
    RobustHelix const& rhel = hseed._helix;
    ComboHitCollection const& hhits = hseed._hhits;

    static unsigned igraph = 0;
    igraph++;
    art::ServiceHandle<art::TFileService> tfs;
    char pri_used_name[100];
    snprintf(pri_used_name,100,"pri_used_shfz%i",igraph);
    char pri_notused_name[100];
    snprintf(pri_notused_name,100,"pri_notused_shfz%i",igraph);
    char bkg_used_name[100];
    snprintf(bkg_used_name,100,"bkg_used_shfz%i",igraph);
    char notselected_name[100];
    snprintf(notselected_name,100,"notselected_shfz%i",igraph);
    char selected_name[100];
    snprintf(selected_name,100,"selected_shfz%i",igraph);
    char tc_name[100];
    snprintf(tc_name,100,"tc_shfz%i",igraph);
    char trk_name[100];
    snprintf(trk_name,100,"trk_shfz%i",igraph);
    char title[100];
    snprintf(title,100,"StrawHit #PhiZ evt %i hel %i;mm;rad",_iev,ihel);
    TH2F* pri_used = tfs->make<TH2F>(pri_used_name,title,100,-_zsize,_zsize,100,-_fsize,_fsize);
    TH2F* pri_notused = tfs->make<TH2F>(pri_notused_name,title,100,-_zsize,_zsize,100,-_fsize,_fsize);
    TH2F* bkg_used = tfs->make<TH2F>(bkg_used_name,title,100,-_zsize,_zsize,100,-_fsize,_fsize);
    TH2F* notselected = tfs->make<TH2F>(notselected_name,title,100,-_zsize,_zsize,100,-_fsize,_fsize);
    TH2F* selected = tfs->make<TH2F>(selected_name,title,100,-_zsize,_zsize,100,-_fsize,_fsize);
    TH2F* tc = tfs->make<TH2F>(tc_name,title,100,-_zsize,_zsize,100,-_fsize,_fsize);
    TH2F* trk = tfs->make<TH2F>(trk_name,title,100,-_zsize,_zsize,100,-_fsize,_fsize);
    pri_used->SetStats(0);
    pri_notused->SetStats(0);
    bkg_used->SetStats(0);
    notselected->SetStats(0);
    selected->SetStats(0);
    tc->SetStats(0);
    trk->SetStats(0);
// get the 'function' lists for these (can be any TObject)
    TList* pri_used_list = pri_used->GetListOfFunctions();
    TList* pri_notused_list = pri_notused->GetListOfFunctions();
    TList* bkg_used_list = bkg_used->GetListOfFunctions();
    TList* notselected_list = notselected->GetListOfFunctions();
    TList* selected_list = selected->GetListOfFunctions();
    TList* tc_list = tc->GetListOfFunctions();
//    TList* trk_list = trk->GetListOfFunctions();
    // MC true hit positions
    TH2F* mct(0);
    TList* mct_list(0);
    if(_mcdiag){
      char mctruth_name[100];
      snprintf(mctruth_name,100,"mctshfz%i",igraph);
      mct = tfs->make<TH2F>(mctruth_name,title,100,-_zsize,_zsize,100,-_fsize,_fsize);
      mct_list = mct->GetListOfFunctions();
    }

    for(size_t ihit=0;ihit < hhits.size(); ++ihit) {
      ComboHit const& hhit = hhits[ihit];
    // compute the projected phi (y) error for this elipse.  Z error is always the straw width
    // create an elipse for this hit
      XYZVec that = XYZVec(-(hhit.pos().y()-rhel.centery()),(hhit.pos().x()-rhel.centerx()),0.0).unit(); // local phi direction at this hit
      double proj = that.Dot(hhit.wdir());
      double fcent = rhel.circleAzimuth(hhit.pos().z());
      TEllipse* te = new TEllipse(hhit.pos().z(),hhit.helixPhi()-fcent,
	hhit._tres, hhit._wres*proj/rhel.radius());

      // mc truth
      if(_mcdiag){
	std::vector<StrawDigiIndex> sdis;
	hhits.fillStrawDigiIndices(evt,ihit,sdis);

	if(primary(pspp,sdis[0])){

	  if (use(hhit) ) {
	    te->SetFillColor(kRed);
	    te->SetLineColor(kRed);
	    pri_used_list->Add(te);
	  } else {
	    te->SetFillColor(kCyan);
	    te->SetLineColor(kCyan);
	    pri_notused_list->Add(te);
	  }
	}
	else {
	  if (use(hhit)) {
	    te->SetFillColor(kMagenta);
	    te->SetLineColor(kMagenta);
	    bkg_used_list->Add(te);
	  } 
	}
      }
    }
  // TC hits
    auto timec = hseed.timeCluster();
    if(timec.isNonnull()){
      for(auto ih : timec->hits()){
	auto const& ch = _chcol->at(ih);
	XYZVec rpos(ch.pos().x()-rhel.centerx(), ch.pos().y()-rhel.centery(), ch.pos().z());
	// resolve phi
	double dphi = atan2(rpos.y(),rpos.x()) - rhel.circleAzimuth(rpos.z());
	dphi -= rint(dphi/CLHEP::twopi)*CLHEP::twopi;
	XYZVec that = XYZVec(-rpos.y(),rpos.x(),0.0).unit(); // local phi direction at this hit
	double proj = that.Dot(ch.wdir());
	TEllipse* te = new TEllipse(ch.pos().z(),dphi,
	    ch._tres, ch._wres*proj/rhel.radius());
	te->SetLineColor(kGreen);
	te->SetFillColor(kGreen);
	tc_list->Add(te);
      }
    }
  // background hits
    for(size_t ich = 0; ich< _chcol->size();++ich) {
      auto const& ch = _chcol->at(ich);
      XYZVec rpos(ch.pos().x()-rhel.centerx(), ch.pos().y()-rhel.centery(), ch.pos().z());
      // resolve phi
      double dphi = atan2(rpos.y(),rpos.x()) - rhel.circleAzimuth(rpos.z());
      dphi -= rint(dphi/CLHEP::twopi)*CLHEP::twopi;
      XYZVec that = XYZVec(-rpos.y(),rpos.x(),0.0).unit(); // local phi direction at this hit
      double proj = that.Dot(ch.wdir());
      TEllipse* te = new TEllipse(ch.pos().z(),dphi,
	  ch._tres, ch._wres*proj/rhel.radius());
      te->SetLineColor(kYellow);
      te->SetFillColor(kYellow);
      notselected_list->Add(te);
      if(selectedHit(ich)){
	TEllipse* tes = new TEllipse(ch.pos().z(),dphi,
	  ch._tres, ch._wres*proj/rhel.radius());
	tes->SetLineColor(kOrange);
	tes->SetFillColor(kOrange);
	selected_list->Add(tes);
      }
      if(_mcdiag){
	// get digi pointers from combo hits
	std::vector<StrawDigiIndex> sdis;
	_chcol->fillStrawDigiIndices(evt,ich,sdis);
	for(auto isd : sdis) { 
	  StrawDigiMC const& mcdigi = _mcdigis->at(isd);
	  art::Ptr<StepPointMC> spmcp;
	  if (TrkMCTools::stepPoint(spmcp,mcdigi) >= 0 && spmcp->simParticle() == pspp){
	    double mcdphi = atan2(spmcp->position().y()-_mch.centery(),spmcp->position().x()-_mch.centerx()) - _mch.circleAzimuth(spmcp->position().z());
	    mcdphi -= rint(mcdphi/CLHEP::twopi)*CLHEP::twopi;
	    TMarker* mch = new TMarker(spmcp->position().z(),mcdphi,20);
	    mch->SetMarkerColor(kBlue);
	    mch->SetMarkerSize(1);
	    mct_list->Add(mch);
	    if(!selectedHit(ich)){
	      te->SetLineColor(kBlue);
	      te->SetFillColor(kBlue);
	    }
	  }
	}
      }
    }
  }

  void HelixDiag::plotXY(art::Event const& evt, art::Ptr<SimParticle> const& pspp, HelixSeed const& hseed, unsigned ihel) {
    RobustHelix const& rhel = hseed._helix;
    ComboHitCollection const& hhits = hseed._hhits;

    static std::string pri_used_title("Pri Used;x(mm);y(mm)");
    static std::string pri_notused_title("Pri Not Used;x(mm);y(mm)");
    static std::string bkg_used_title("Bkg Used;x(mm);y(mm)");
    static std::string notselected_title("Not Selected;x(mm);y(mm)");
    static std::string selected_title("Selected;x(mm);y(mm)");
    static std::string tc_title("Time Cluster;x(mm);y(mm)");
    static std::string trk_title("Tracker XY;x(mm);y(mm)");
    static std::string mctruth_title("MC Truth;x(mm);y(mm)");

    static unsigned igraph = 0;
    igraph++;
    art::ServiceHandle<art::TFileService> tfs;
    char pri_used_name[100];
    snprintf(pri_used_name,100,"pri_used_shxy%i",igraph);
    char pri_notused_name[100];
    snprintf(pri_notused_name,100,"pri_notused_shxy%i",igraph);
    char bkg_used_name[100];
    snprintf(bkg_used_name,100,"bkg_used_shxy%i",igraph);
    char notselected_name[100];
    snprintf(notselected_name,100,"notselected_shxy%i",igraph);
    char selected_name[100];
    snprintf(selected_name,100,"selected_shxy%i",igraph);
    char tc_name[100];
    snprintf(tc_name,100,"tc_shxy%i",igraph);
    char trk_name[100];
    snprintf(trk_name,100,"trk_shxy%i",igraph);
    char title[100];
    snprintf(title,100,"StrawHit XY evt %i hel %i;mm;rad",_iev,ihel);
    TH2F* pri_used = tfs->make<TH2F>(pri_used_name,pri_used_title.c_str(),100,-_xysize,_xysize,100,-_xysize,_xysize);
    TH2F* pri_notused = tfs->make<TH2F>(pri_notused_name,pri_notused_title.c_str(),100,-_xysize,_xysize,100,-_xysize,_xysize);
    TH2F* bkg_used = tfs->make<TH2F>(bkg_used_name,bkg_used_title.c_str(),100,-_xysize,_xysize,100,-_xysize,_xysize);
    TH2F* notselected = tfs->make<TH2F>(notselected_name,notselected_title.c_str(),100,-_xysize,_xysize,100,-_xysize,_xysize);
    TH2F* selected = tfs->make<TH2F>(selected_name,selected_title.c_str(),100,-_xysize,_xysize,100,-_xysize,_xysize);
    TH2F* tc = tfs->make<TH2F>(tc_name,tc_title.c_str(),100,-_xysize,_xysize,100,-_xysize,_xysize);
    TH2F* trk = tfs->make<TH2F>(trk_name,trk_title.c_str(),100,-_xysize,_xysize,100,-_xysize,_xysize);
    pri_used->SetStats(0);
    pri_notused->SetStats(0);
    bkg_used->SetStats(0);
    notselected->SetStats(0);
    selected->SetStats(0);
    tc->SetStats(0);
    trk->SetStats(0);

    pri_used->SetLineColor(kRed);
    pri_notused->SetLineColor(kCyan);
    bkg_used->SetLineColor(kMagenta);
    notselected->SetLineColor(kYellow);
    selected->SetLineColor(kOrange);
    tc->SetLineColor(kGreen);
// get the 'function' lists for these (can be any TObject)
    TList* pri_used_list = pri_used->GetListOfFunctions();
    TList* pri_notused_list = pri_notused->GetListOfFunctions();
    TList* bkg_used_list = bkg_used->GetListOfFunctions();
    TList* notselected_list = notselected->GetListOfFunctions();
    TList* selected_list = selected->GetListOfFunctions();
    TList* tc_list = tc->GetListOfFunctions();
    TList* trk_list = trk->GetListOfFunctions();
    // MC true hit positions
    TH2F* mct(0);
    TList* mct_list(0);
    if(_mcdiag){
      char mctruth_name[100];
      snprintf(mctruth_name,100,"mctshxy%i",igraph);
      mct = tfs->make<TH2F>(mctruth_name,mctruth_title.c_str(),100,-_xysize,_xysize,100,-_xysize,_xysize);
      mct->SetStats(0);
      mct->SetLineColor(kBlue);
      mct->SetMarkerStyle(20);
      mct->SetMarkerColor(kBlue);
      mct_list = mct->GetListOfFunctions();
    }
    // define plot center
    XYZVec pcent;
    if(_hcenter) pcent = XYZVec( rhel.centerx(), rhel.centery(), 0.0);
    // draw the detector boundaries
    static double innerrad(380.0);
    static double outerrad(680.0);
    TArc* indet = new TArc(-pcent.x(),-pcent.y(),innerrad);
    TArc* outdet = new TArc(-pcent.x(),-pcent.y(),outerrad);
    indet->SetLineColor(kBlack);
    indet->SetFillColor(kWhite);
    indet->SetFillStyle(1001);
    outdet->SetLineColor(kBlack);
    outdet->SetFillColor(kGray);
//    outdet->SetFillStyle(3444);

    TArc* target = new TArc(-pcent.x(),-pcent.y(),_targetradius);
    target->SetLineColor(kBlack);
    target->SetFillColor(kBlack);
    target->SetFillStyle(3002);
    // add these to the plot
    trk_list->Add(outdet);
    trk_list->Add(indet);
    trk_list->Add(target);

    for(size_t ihit=0;ihit <  hhits.size(); ++ihit ) {
      ComboHit const& hhit = hhits[ihit];
    // create an elipse for this hit
      TEllipse* te = new TEllipse(hhit.pos().x()-pcent.x(),hhit.pos().y()-pcent.y(),
	  hhit._wres, hhit._tres,0.0,360.0, hhit.wdir().phi()*180.0/TMath::Pi());
      // mc truth
      if(_mcdiag){
	std::vector<StrawDigiIndex> sdis;
	hhits.fillStrawDigiIndices(evt,ihit,sdis);

	if(primary(pspp,sdis[0])){
	  if (use(hhit) ) {
	    te->SetFillColor(kRed);
	    te->SetLineColor(kRed);
	    pri_used_list->Add(te);
	  } else {
	    te->SetFillColor(kCyan);
	    te->SetLineColor(kCyan);
	    pri_notused_list->Add(te);
	  }
	} else {
	  if (use(hhit)) {
	    te->SetFillColor(kMagenta);
	    te->SetLineColor(kMagenta);
	    bkg_used_list->Add(te);
	  } 
	}
      }
    }
  // TC hits
    auto timec = hseed.timeCluster();
    if(timec.isNonnull()){
      for(auto ih : timec->hits()){
	auto const& ch = _chcol->at(ih);
	TEllipse* te = new TEllipse(ch._pos.x()-pcent.x(),ch._pos.y()-pcent.y(),
	    ch._wres, ch._tres,0.0,360.0, ch._wdir.phi()*180.0/TMath::Pi());
	te->SetLineColor(kGreen);
	te->SetFillColor(kGreen);
	tc_list->Add(te);
      }
    }
  // all hits
    for(size_t ich = 0; ich< _chcol->size();++ich) {
      auto const& ch = _chcol->at(ich);
      // create an elipse for this hit
      TEllipse* te = new TEllipse(ch._pos.x()-pcent.x(),ch._pos.y()-pcent.y(),
	  ch._wres, ch._tres,0.0,360.0, ch._wdir.phi()*180.0/TMath::Pi());
      // draw all hits
      te->SetLineColor(kYellow);
      te->SetFillColor(kYellow);
      notselected_list->Add(te);
      // selected hits; this overlaps with those below
      if(selectedHit(ich)){
	TEllipse* tes = new TEllipse(ch._pos.x()-pcent.x(),ch._pos.y()-pcent.y(),
	  ch._wres, ch._tres,0.0,360.0, ch._wdir.phi()*180.0/TMath::Pi());
 	tes->SetLineColor(kOrange);
	tes->SetFillColor(kOrange);
	selected_list->Add(tes);
      }
      if(_mcdiag){
	// get digi pointers from combo hits
	std::vector<StrawDigiIndex> sdis;
	_chcol->fillStrawDigiIndices(evt,ich,sdis);
	for(auto isd : sdis) { 
	  StrawDigiMC const& mcdigi = _mcdigis->at(isd);
	  art::Ptr<StepPointMC> spmcp;
	  if (TrkMCTools::stepPoint(spmcp,mcdigi) >= 0 && spmcp->simParticle() == pspp){
	    TMarker* mch = new TMarker(spmcp->position().x()-pcent.x(),spmcp->position().y()-pcent.y(),20);
	    mch->SetMarkerColor(kBlue);
	    mch->SetMarkerSize(1);
	    mct_list->Add(mch);
	    if(!selectedHit(ich)){
	      te->SetLineColor(kBlue);
	      te->SetFillColor(kBlue);
	    }
	  }
	}
      }
    }
    // add fit circle
    TArc* fitarc = new TArc(rhel.centerx()-pcent.x(),rhel.centery()-pcent.y(),rhel.radius());
    fitarc->SetLineColor(kRed);
    fitarc->SetLineWidth(2);
    fitarc->SetFillStyle(0);
    pri_used_list->Add(fitarc);
    if (_mcdiag) {
      // mc true circle
      TArc* mcarc = new TArc(_mch.centerx()-pcent.x(),_mch.centery()-pcent.y(),_mch.radius());
      mcarc->SetLineColor(kBlue);
      mcarc->SetLineWidth(2);
      mcarc->SetFillStyle(0);
      mct_list->Add(mcarc);
    }
  }

   void HelixDiag::fillHitInfoMC(const art::Ptr<SimParticle>& pspp, StrawDigiMC const& digimc, HitInfoMC& hinfomc) {
    hinfomc.reset();
    art::Ptr<SimParticle> spp;
    art::Ptr<StepPointMC> spmcp;
    if(TrkMCTools::simParticle(spp,digimc) > 0 && TrkMCTools::stepPoint(spmcp,digimc) >= 0 ){
      hinfomc._pdg = spp->pdgId();
      hinfomc._proc = spp->originParticle().creationCode();
      if(spp->genParticle().isNonnull())
	hinfomc._gen = spp->genParticle()->generatorId().id();
      MCRelationship rel(pspp,spp);
      hinfomc._rel = rel.relationship();
      hinfomc._t0 = _toff.timeWithOffsetsApplied(*spmcp);
    }
  }

  bool HelixDiag::fillMCHelix(art::Ptr<SimParticle> const& pspp) {
    bool retval(false);
    GeomHandle<DetectorSystem> det;
    GlobalConstantsHandle<ParticleDataTable> pdt;
    double charge = pdt->particle(pspp->pdgId()).ref().charge();
    Hep3Vector pos;
    Hep3Vector mom;
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
      mom = jmc->momentum();
      pos = det->toDetector(jmc->position());
      retval = true;
    } else {
    // no vd steps: try to fill from the tracker StepPoints
      for(auto imcd = _mcdigis->begin();imcd != _mcdigis->end(); ++imcd){
	auto const& stepptr = imcd->stepPointMC(imcd->earlyEnd());
	if(stepptr->simParticle() == pspp){
	  mom = stepptr->momentum();
	  pos = stepptr->position();
	  retval = true;
	  break;
	}
      }
    }
    if(retval){
      _mcmom = mom.mag();
      _mcpz = mom.z();
      TrkUtilities::RobustHelixFromMom(pos,mom,charge,_bz0,_mch);
    } else {
      _mcmom = -1.0;
      _mcpz = 0.0;
      _mch = RobustHelix();
    }
    return retval;
  }

  bool
  HelixDiag::use(ComboHit const& hhit) const {
     return !hhit._flag.hasAnyProperty(_dontuseflag);
  }

  bool
  HelixDiag::stereo(ComboHit const& hhit) const {
    static StrawHitFlag stereo(StrawHitFlag::stereo);
    return hhit._flag.hasAllProperties(stereo);
  }

  bool
  HelixDiag::primary(art::Ptr<SimParticle> const& pspp, size_t index) {
    StrawDigiMC const& mcdigi = _mcdigis->at(index);
    art::Ptr<StepPointMC> spmcp;
    return TrkMCTools::stepPoint(spmcp,mcdigi) >= 0 && spmcp->simParticle() == pspp;
  }

  bool
  HelixDiag::selectedHit(size_t index) {
    StrawHitFlag flag;
    if(_useshfcol)
      flag = _shfcol->at(index);
    else
      flag = _chcol->at(index).flag();
    return flag.hasAllProperties(_hsel) && !flag.hasAnyProperty(_hbkg);
  }

}

using mu2e::HelixDiag;
DEFINE_ART_MODULE(HelixDiag);
