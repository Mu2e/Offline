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
// diagnostics
#include "TrkDiag/inc/TrkMCTools.hh"
#include "TrkReco/inc/XYZP.hh"
#include "TrkDiag/inc/HitInfoMC.hh"
#include "TrkDiag/inc/TrkMCTools.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/HelixSeedCollection.hh"
#include "DataProducts/inc/threevec.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
// root
#include "TGraph.h"
#include "TH1F.h"
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
      virtual void beginJob();
      // This is called for each event.
      virtual void analyze(art::Event const& e);
    private:
// config parameters
      int _diag;
      bool _mcdiag;
      unsigned _minnce; // minimum # CE hits to make plots
      double _targetradius;
      bool _plotxy, _plotz;
      // event object tags      art::InputTag _shTag;
      art::InputTag _shTag;
      art::InputTag _shpTag;
      art::InputTag _shfTag;
      art::InputTag _hsTag;
      art::InputTag _mcdigisTag;
 
      // cache of event objects
      const StrawHitCollection* _shcol;
      const StrawHitPositionCollection* _shpcol;
      const StrawHitFlagCollection* _shfcol;
      const HelixSeedCollection* _hscol;
      const StrawDigiMCCollection* _mcdigis;
      // helper functions
      bool findData(const art::Event& e);
      bool conversion(size_t index);
      unsigned countCEHits(vector<hitIndex> const& hits );
      void setPhi(RobustHelix const& helix,	XYZPVector& xyzp);
      void fillHitInfoMC(const art::Ptr<SimParticle>& pspp, StrawDigiMC const& digimc, HitInfoMC& hinfomc);
      // display functions
      void plotXY(XYZPVector const& xyzp, HelixSeed const& myseed);
      void plotZ(XYZPVector const& xyzp, HelixSeed const& myseed);
      // TTree and branch variables
      TTree *_hdiag;
      RobustHelix _recoh;
      Int_t _nhits, _nprimary;
      std::vector<HitInfoMC> _hinfomc;
  };

  HelixDiag::~HelixDiag() {
  }

  HelixDiag::HelixDiag(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _diag(pset.get<int>("DiagLevel",1)),
    _mcdiag(pset.get<bool>("MonteCarloDiag",true)),
    _minnce(pset.get<unsigned>("MinimumCEHits",10)),
    _targetradius(pset.get<double>("TargetRadius",75)),
    _plotxy(pset.get<bool>("PlotXY",true)),
    _plotz(pset.get<bool>("PlotZ",true)),
    _shTag(pset.get<string>("StrawHitCollectionTag","makeSH")),
    _shpTag(pset.get<string>("StrawHitPositionCollectionTag","MakeStereoHits")),
    _shfTag(pset.get<string>("StrawHitFlagCollectionTag","PosHelixFinder")),
    _hsTag(pset.get<string>("HelixSeedCollectionTag","PosHelixFinder")),
    _mcdigisTag(pset.get<art::InputTag>("StrawDigiMCCollection","makeSH"))
   {
    if(_diag > 0){
      art::ServiceHandle<art::TFileService> tfs;
      _hdiag=tfs->make<TTree>("hdiag","Helix Finding diagnostics");
      _hdiag->Branch("recoh",&_recoh);
      _hdiag->Branch("nhits",&_nhits,"nhits/I");
      if(_mcdiag){
	_hdiag->Branch("nprimary",&_nprimary,"nprimary/I");
	_hdiag->Branch("tshmc",&_hinfomc);
      }
    }
  }

  void HelixDiag::beginJob(){
  }

  void HelixDiag::analyze(art::Event const& evt) {
// find the data
    if(findData(evt)) {
    // loop over helices
      for(auto hseed : *_hscol) {
	RobustHelix const& myhel = hseed._helix;
	vector<hitIndex> const& hits = hseed._timeCluster._strawHitIdxs; 
	// fill XYZP points used in this helix
	XYZPVector xyzp;
	XYZP::fillXYZP(*_shcol, *_shpcol, hits, xyzp);
	// resolve the phi for these points.  Note that this isn't necessarily the same resolution
	// as used in the original fit
	setPhi(myhel,xyzp);
	// fill TTree branches
	_recoh = myhel;
	_nhits = hits.size();
	art::Ptr<SimParticle> pspp;
	TrkMCTools::primaryParticle(pspp,hits,_mcdigis);
	for( auto hi : hits ) {
	  StrawDigiMC const& digimc = _mcdigis->at(hi);
	  HitInfoMC hinfomc;
	  fillHitInfoMC(pspp,digimc,hinfomc);
	  _hinfomc.push_back(hinfomc);
	  
	  
	}
	// fill the tree
	_hdiag->Fill();

	if(_plotxy || _plotz ) {
	// cou(nt the # of conversion hits in this helix
	  unsigned nce = countCEHits(hits);
	  if (nce >= _minnce) {
	    // fill graphs for display
	    if(_plotxy)plotXY(xyzp,hseed);
	    if(_plotz)plotZ(xyzp,hseed);
	  }
	}
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
    }

    return _shcol != 0 && _shpcol != 0 && _shfcol != 0 && _hscol != 0 && (_mcdigis != 0 || !_mcdiag);
  }


  void HelixDiag::plotXY(XYZPVector const& xyzp, HelixSeed const& hseed) {
    RobustHelix const& myhel = hseed._helix;
    vector<hitIndex> const& hits = hseed._timeCluster._strawHitIdxs; 

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
    snprintf(title,100,"StrawHit XY trk %i;mm;rad",igraph);
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

    for(unsigned ih=0;ih<xyzp.size();++ih){
      if(conversion(static_cast<size_t>(xyzp[ih]._ind))){
	if (xyzp[ih].use()) {
	  if (xyzp[ih].stereo()) {
	    ce_stereo_used->Fill(xyzp[ih]._pos.x()-myhel.center().x(),xyzp[ih]._pos.y()-myhel.center().y());
	  }
	  else {
	    ce_notstereo_used->Fill(xyzp[ih]._pos.x()-myhel.center().x(),xyzp[ih]._pos.y()-myhel.center().y());
	  }	      
	}
	else {
	  if (xyzp[ih].stereo()) {
	    ce_stereo_notused->Fill(xyzp[ih]._pos.x()-myhel.center().x(),xyzp[ih]._pos.y()-myhel.center().y());
	  }
	  else {
	    ce_notstereo_notused->Fill(xyzp[ih]._pos.x()-myhel.center().x(),xyzp[ih]._pos.y()-myhel.center().y());
	  }
	}
      }
      else {
	if (xyzp[ih].use()) {
	  if (xyzp[ih].stereo()) {
	    bkg_stereo_used->Fill(xyzp[ih]._pos.x()-myhel.center().x(),xyzp[ih]._pos.y()-myhel.center().y());
	  }
	  else {
	    bkg_notstereo_used->Fill(xyzp[ih]._pos.x()-myhel.center().x(),xyzp[ih]._pos.y()-myhel.center().y());
	  }	      
	}
	else {
	  if (xyzp[ih].stereo()) {
	    bkg_stereo_notused->Fill(xyzp[ih]._pos.x()-myhel.center().x(),xyzp[ih]._pos.y()-myhel.center().y());
	  }
	  else {
	    bkg_notstereo_notused->Fill(xyzp[ih]._pos.x()-myhel.center().x(),xyzp[ih]._pos.y()-myhel.center().y());
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

      for(auto hit : hits ) {
	StrawDigiMC const& mcdigi = _mcdigis->at(hit);
	art::Ptr<StepPointMC> spmcp;
	if (TrkMCTools::stepPoint(spmcp,mcdigi) >= 0 && TrkMCTools::CEDigi(mcdigi)){
	  mct->Fill(spmcp->position().x()-myhel.center().x(),spmcp->position().y()-myhel.center().y());
	}
      }
    }
  }

  void HelixDiag::plotZ(XYZPVector const& xyzp, HelixSeed const& hseed) {
    RobustHelix const& myhel = hseed._helix;
    vector<hitIndex> const& hits = hseed._timeCluster._strawHitIdxs; 
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
    snprintf(title,100,"StrawHit #phi Z trk %i;mm;rad",igraph);
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

    for(unsigned ih=0;ih<xyzp.size();++ih){
      if(conversion(static_cast<size_t>(xyzp[ih]._ind))){
	if (xyzp[ih].use()) {
	  if (xyzp[ih].stereo()) {
	    ce_stereo_used->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	  }
	  else {
	    ce_notstereo_used->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	  }	      
	}
	else {
	  if (xyzp[ih].stereo()) {
	    ce_stereo_notused->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	  }
	  else {
	    ce_notstereo_notused->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	  }
	}
      }
      else {
	if (xyzp[ih].use()) {
	  if (xyzp[ih].stereo()) {
	    bkg_stereo_used->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	  }
	  else {
	    bkg_notstereo_used->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	  }	      
	}
	else {
	  if (xyzp[ih].stereo()) {
	    bkg_stereo_notused->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
	  }
	  else {
	    bkg_notstereo_notused->Fill(xyzp[ih]._pos.z(),xyzp[ih]._phi);
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

      for(auto& hit : hits ) {
	StrawDigiMC const& mcdigi = _mcdigis->at(hit);
	art::Ptr<StepPointMC> spmcp;
	if (TrkMCTools::stepPoint(spmcp,mcdigi) >= 0 && TrkMCTools::CEDigi(mcdigi)){
	  mct->Fill(spmcp->position().z(),spmcp->position().phi());
	}
      }
    }
  }

  unsigned HelixDiag::countCEHits(vector<hitIndex> const& hitindexs ) {
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

  void HelixDiag::setPhi(RobustHelix const& helix, XYZPVector& xyzpv){
// compare expected phi position with actual, and adjust the phase to make these consistent
    for(auto& xyzp : xyzpv) {
      // compute the expected phi from the z position of the straw
      double phiex = helix.fz0() + xyzp._pos.z()/helix.lambda();
      // compute phi WRT the circle center from the transverse position
      double phi = Hep3Vector(xyzp._pos - helix.center()).phi();
      // minimize the phase within the 2pi ambiguity
      double dphi = Angles::deltaPhi(phi,phiex);
      if(fabs(dphi) > M_PI ) std::cout << "Failed to resolve phi" << std::endl;
      xyzp._phi = phi;
    }

  }

  void HelixDiag::fillHitInfoMC(const art::Ptr<SimParticle>& pspp, StrawDigiMC const& digimc, HitInfoMC& hinfomc) {
    art::Ptr<SimParticle> spp;
    hinfomc.reset();
    if(TrkMCTools::simParticle(spp,digimc) > 0){
      hinfomc._pdg = spp->pdgId();
      hinfomc._proc = spp->realCreationCode();
      if(spp->genParticle().isNonnull())
	hinfomc._gen = spp->genParticle()->generatorId().id();
      hinfomc._rel = MCRelationship::relationship(pspp,spp);
    }
  }

}

using mu2e::HelixDiag;
DEFINE_ART_MODULE(HelixDiag);
