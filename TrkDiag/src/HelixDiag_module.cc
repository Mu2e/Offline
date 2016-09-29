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
#include "TrkDiag/inc/HitInfoMC.hh"
#include "TrkDiag/inc/TrkMCTools.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/HelixSeedCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
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
      unsigned _minnce; // minimum # CE hits to make plots
      double _targetradius;
      bool _plot;
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
      // time offsets
      SimParticleTimeOffset _toff;
      // Virtual Detector IDs
      std::vector<int> _midvids;
      // cache of BField at 0,0,0
      double _bz0;
      // helper functions
      bool findData(const art::Event& e);
      bool conversion(size_t index);
      unsigned countCEHits(vector<StrawHitIndex> const& hits );
      void fillHitInfoMC(art::Ptr<SimParticle> const& pspp, StrawDigiMC const& digimc, HitInfoMC& hinfomc);
      void fillMCHelix(art::Ptr<SimParticle> const& pspp);
      // display functions
      void plotXY(HelixSeed const& myseed, unsigned ihel);
      void plotZ(HelixSeed const& myseed, unsigned ihel);
      bool use(HelixHit const& hhit) const;
      bool stereo(HelixHit const& hhit) const;
       // TTree and branch variables
      TTree *_hdiag;
      Int_t _iev;
      RobustHelix _reh;
      Int_t _nhits, _nprimary;
      Int_t _pdg, _gen, _proc;
      std::vector<HitInfoMC> _hinfomc;
      RobustHelix _mch;
      Hep3Vector _mcmom;
  };

  HelixDiag::~HelixDiag() {
  }

  HelixDiag::HelixDiag(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _diag(pset.get<int>("DiagLevel",1)),
    _dontuseflag(pset.get<std::vector<std::string>>("UseFlag",vector<string>{"Outlier","OtherBackground"})),
    _mcdiag(pset.get<bool>("MonteCarloDiag",true)),
    _minnce(pset.get<unsigned>("MinimumCEHits",10)),
    _targetradius(pset.get<double>("TargetRadius",75)),
    _plot(pset.get<bool>("PlotHelices",false)),
    _shTag(pset.get<string>("StrawHitCollectionTag","makeSH")),
    _shpTag(pset.get<string>("StrawHitPositionCollectionTag","MakeStereoHits")),
    _shfTag(pset.get<string>("StrawHitFlagCollectionTag","PosHelixFinder")),
    _hsTag(pset.get<string>("HelixSeedCollectionTag","PosHelixFinder")),
    _mcdigisTag(pset.get<art::InputTag>("StrawDigiMCCollection","makeSH")),
    _vdmcstepsTag(pset.get<art::InputTag>("VDStepPointMCCollection","detectorFilter:virtualdetector")),
    _toff(pset.get<fhicl::ParameterSet>("TimeOffsets"))
  {
    if(_diag > 0){
      art::ServiceHandle<art::TFileService> tfs;
      _hdiag=tfs->make<TTree>("hdiag","Helix Finding diagnostics");
      _hdiag->Branch("iev",&_iev,"iev/I");
      _hdiag->Branch("reh",&_reh);
      _hdiag->Branch("nhits",&_nhits,"nhits/I");
      if(_mcdiag){
	_hdiag->Branch("mch",&_mch);
	_hdiag->Branch("mcmom",&_mcmom);
	_hdiag->Branch("nprimary",&_nprimary,"nprimary/I");
	_hdiag->Branch("pdg",&_pdg,"pdg/I");
	_hdiag->Branch("gen",&_gen,"gen/I");
	_hdiag->Branch("proc",&_proc,"proc/I");
	if(_diag > 1){
	  _hdiag->Branch("tshmc",&_hinfomc);
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
      for(auto hseed : *_hscol) {
	RobustHelix const& myhel = hseed._helix;
	HelixHitCollection& hhits = hseed._hhits;
	std::vector<StrawHitIndex> hits;
	for(auto hhit : hhits) {
	  hits.push_back(hhit._shidx);
	}
	// fill TTree branches; first the ones that don't depend on MC
	_reh = myhel;
	if(_mcdiag) {
	  // get information about the primary particle (produced most hits)
	  art::Ptr<SimParticle> pspp;
	  _nprimary = TrkMCTools::primaryParticle(pspp,hits,_mcdigis);
	  _pdg = pspp->pdgId();
	  _proc = pspp->realCreationCode();
	  if(pspp->genParticle().isNonnull())
	    _gen = pspp->genParticle()->generatorId().id();
	  else
	    _gen = -1;
	  // fill MC true helix parameters
	  fillMCHelix(pspp);
	  if(_diag > 1){
	// fill information about individual hits
	    for( auto hi : hits ) {
	      StrawDigiMC const& digimc = _mcdigis->at(hi);
	      HitInfoMC hinfomc;
	      fillHitInfoMC(pspp,digimc,hinfomc);
	      _hinfomc.push_back(hinfomc);
	    }
	  }
	}
	  // fill the tree
	_hdiag->Fill();
	// if requested, plot the hits and helices
	if( _plot ) {
	// count the # of conversion hits in this helix
	  unsigned nce = countCEHits(hits);
	  if (nce >= _minnce) {
	    // fill graphs for display
	    plotXY(hseed,ihel);
	    plotZ(hseed,ihel);
	  }
	}
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


  void HelixDiag::plotXY(HelixSeed const& hseed, unsigned ihel) {
    RobustHelix const& myhel = hseed._helix;
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

    for(auto hhit : hhits ) {
      if(conversion(hhit._shidx)){
	if (use(hhit) ) {
	  if (stereo(hhit)) {
	    ce_stereo_used->Fill(hhit._pos.x()-myhel.center().x(),hhit._pos.y()-myhel.center().y());
	  }
	  else {
	    ce_notstereo_used->Fill(hhit._pos.x()-myhel.center().x(),hhit._pos.y()-myhel.center().y());
	  }	      
	}
	else {
	  if (stereo(hhit)) {
	    ce_stereo_notused->Fill(hhit._pos.x()-myhel.center().x(),hhit._pos.y()-myhel.center().y());
	  }
	  else {
	    ce_notstereo_notused->Fill(hhit._pos.x()-myhel.center().x(),hhit._pos.y()-myhel.center().y());
	  }
	}
      }
      else {
	if (use(hhit)) {
	  if (stereo(hhit)) {
	    bkg_stereo_used->Fill(hhit._pos.x()-myhel.center().x(),hhit._pos.y()-myhel.center().y());
	  }
	  else {
	    bkg_notstereo_used->Fill(hhit._pos.x()-myhel.center().x(),hhit._pos.y()-myhel.center().y());
	  }	      
	}
	else {
	  if (stereo(hhit)) {
	    bkg_stereo_notused->Fill(hhit._pos.x()-myhel.center().x(),hhit._pos.y()-myhel.center().y());
	  }
	  else {
	    bkg_notstereo_notused->Fill(hhit._pos.x()-myhel.center().x(),hhit._pos.y()-myhel.center().y());
	  }
	}
      }
    }


    TArc* fitarc = new TArc(0.0,0.0,myhel.radius());
    fitarc->SetLineColor(kRed);
    fitarc->SetLineWidth(2);
    fitarc->SetFillStyle(0);
    // draw the detector boundaries
    static double innerrad(380.0);
    static double outerrad(680.0);
    TArc* indet = new TArc(-myhel.center().x(),-myhel.center().y(),innerrad);
    TArc* outdet = new TArc(-myhel.center().x(),-myhel.center().y(),outerrad);
    indet->SetLineColor(kBlue);
    indet->SetFillStyle(0);
    outdet->SetLineColor(kBlue);
    outdet->SetFillStyle(0);

    TArc* target = new TArc(-myhel.center().x(),-myhel.center().y(),_targetradius);
    target->SetLineColor(kBlack);
    target->SetFillStyle(0);
    // add these to the plot
    TList* flist = ce_stereo_used->GetListOfFunctions();
    flist->Add(fitarc);
    flist->Add(indet);
    flist->Add(outdet);
    flist->Add(target);

    if (_mcdiag) {
      // Plot the MC true CE hit positions
      char mctruth_name[100];
      snprintf(mctruth_name,100,"mctshxy%i",igraph);
      TH2F* mct = tfs->make<TH2F>(mctruth_name,title,100,-500,500,100,-500,500);
      mct->SetMarkerStyle(5);
      mct->SetMarkerColor(kMagenta);

      for(auto hhit : hhits ) {
	StrawDigiMC const& mcdigi = _mcdigis->at(hhit._shidx);
	art::Ptr<StepPointMC> spmcp;
	if (TrkMCTools::stepPoint(spmcp,mcdigi) >= 0 && TrkMCTools::CEDigi(mcdigi)){
	  mct->Fill(spmcp->position().x()-myhel.center().x(),spmcp->position().y()-myhel.center().y());
	}
      }
    }
  }

  void HelixDiag::plotZ(HelixSeed const& hseed, unsigned ihel) {
    RobustHelix const& myhel = hseed._helix;
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

    for(auto hhit : hhits) {
      if(conversion(hhit._shidx)){
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
    TF1* line = new TF1("line","[0]+[1]*x",-1500,1500);
    line->SetParameter(0,myhel.fz0());
    line->SetParameter(1,1.0/myhel.lambda());
    line->SetLineColor(kRed);
    TList* flist = ce_stereo_used->GetListOfFunctions();
    flist->Add(line);

    if (_mcdiag) {
      // Plot the MC true CE hits
      char mctruth_name[100];
      snprintf(mctruth_name,100,"mctshphiz%i",igraph);
      TH2F* mct = tfs->make<TH2F>(mctruth_name,title,100,-1500,1500,100,-12.5,12.5);
      mct->SetMarkerStyle(5);
      mct->SetMarkerColor(kMagenta);

      for(auto& hhit : hhits ) {
	StrawDigiMC const& mcdigi = _mcdigis->at(hhit._shidx);
	art::Ptr<StepPointMC> spmcp;
	if (TrkMCTools::stepPoint(spmcp,mcdigi) >= 0 && TrkMCTools::CEDigi(mcdigi)){
	  mct->Fill(spmcp->position().z(),spmcp->position().phi());
	}
      }
    }
  }

  unsigned HelixDiag::countCEHits(vector<StrawHitIndex> const& hitindexs ) {
    unsigned retval(0);
    for(auto hi : hitindexs ) {
      if(conversion(static_cast<size_t>(hi)))++retval; 
    }
    return retval;
  }

  bool HelixDiag::conversion(size_t index) {
    StrawDigiMC const& mcdigi = _mcdigis->at(index);
    return TrkMCTools::CEDigi(mcdigi);
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

  void HelixDiag::fillMCHelix(art::Ptr<SimParticle> const& pspp) {
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
      _mcmom = jmc->momentum();
      Hep3Vector pos = det->toDetector(jmc->position());
      double charge = pdt->particle(pspp->pdgId()).ref().charge();
      TrkUtilities::RobustHelixFromMom(pos,_mcmom,charge,_bz0,_mch);
    }
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


}

using mu2e::HelixDiag;
DEFINE_ART_MODULE(HelixDiag);
