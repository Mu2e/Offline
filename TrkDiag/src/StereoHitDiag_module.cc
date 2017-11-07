//
// Stereo hit diagnostics.  Split out of MakeStereoHits
//
// Original author D. Brown
//

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
// root
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TMarker.h"
#include "TList.h"
#include "TLegend.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TTree.h"
// data
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitPosition.hh"
#include "RecoDataProducts/inc/StereoHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
// Utilities
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "TrkDiag/inc/TrkMCTools.hh"
// diagnostics
#include "TrkDiag/inc/StrawHitInfo.hh"
#include "DataProducts/inc/threevec.hh"
using namespace std; 

namespace mu2e 
{

  class StereoHitDiag : public art::EDAnalyzer {
    public:
      explicit StereoHitDiag(fhicl::ParameterSet const&);
      virtual ~StereoHitDiag();
      virtual void beginJob();
      virtual void analyze(const art::Event& e);
    private:
      // configuration
      int _diag;
      bool _mcdiag;
      bool _drawstations;
      // event object Tags
      art::InputTag   _shTag;
      art::InputTag   _shpTag;
      art::InputTag   _sthpTag;
      art::InputTag   _stTag;
      art::InputTag   _mcdigisTag;
      // event data cache
      const StrawHitCollection* _shcol;
      const StrawHitPositionCollection *_shpcol, *_sthpcol;
      const StereoHitCollection* _stcol;
      const StrawDigiMCCollection *_mcdigis;
      // diagnostics
      TH1F* _nhits;
      TH1F* _deltat;
      TH1F* _deltaE;
      TH1F* _deltaz;
      TH1F* _sep;
      TH1F* _ddoth;
      TH1F* _dperph;
      TH1F* _dL;
      TH1F* _mva;
      vector<TH2F*> _stations;

      TTree *_spdiag, *_sdiag;
      Float_t _shphi, _stphi, _mcshphi;
      Float_t _shrho, _strho, _mcshrho;
      Float_t _de, _dt, _dist, _dperp, _dz, _rho, _dl1, _dl2;
      Bool_t _tdiv1, _tdiv2;
      Float_t _dc1, _dc2, _chi2, _ndof, _mvaout, _ddot;
      Float_t _schi2, _smvaout, _sddot, _sdist, _sdz;
      Float_t _mcdist, _mcperp;
      Int_t _stereo, _fs, _sfs, _mcr, _mcrel, _mcpdg, _mcgen, _mcproc;
      // helper functions
      void drawStations();
      bool findData(const art::Event& evt);
  };

  StereoHitDiag::StereoHitDiag(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _diag		(pset.get<int>("diagLevel",1)),
    _mcdiag		(pset.get<bool>("MCdiag",true)),
    _drawstations     	(pset.get<bool>("DrawStations",false)),
    _shTag		(pset.get<art::InputTag>("StrawHitCollection","makeSH")),
    _shpTag		(pset.get<art::InputTag>("StrawHitPositionCollection","MakeStrawHitPositions")),
    _sthpTag		(pset.get<art::InputTag>("StereoHitPositionCollection","MakeStereoHits")),
    _stTag		(pset.get<art::InputTag>("StereoHitCollection","MakeStereoHits")),
    _mcdigisTag		(pset.get<art::InputTag>("StrawDigiMCCollection","makeSD"))
  {}

  StereoHitDiag::~StereoHitDiag(){}

  void StereoHitDiag::beginJob() {
    // create diagnostics if requested
    if(_diag > 0){
      art::ServiceHandle<art::TFileService> tfs;
      _nhits = tfs->make<TH1F>("nhits","NHits",500,0,5000);
      _deltat = tfs->make<TH1F>("deltat","#Delta t;ns",100,-200.0,200.0);
      _deltaE = tfs->make<TH1F>("deltaE","#Delta E/#Sigma E;Ratio",100,-1.0,1.0);
      _deltaz = tfs->make<TH1F>("deltaz","#Delta d;mm",120,-120.0,120.0);
      _ddoth = tfs->make<TH1F>("ddot","#Delta drection;cos(#theta)",120,-3.15,3.15);
      _dperph = tfs->make<TH1F>("dperp","#Delta #rho;mm",120,0.0,200.0);
      _sep = tfs->make<TH1F>("sep","Face separation",6,-0.5,5.5);
      _dL = tfs->make<TH1F>("dL","Length Difference;mm",100,0.0,700.0);
      _mva = tfs->make<TH1F>("mva","MVA output",100,-0.05,1.05);
      if( _diag > 1){
	// detailed diagnostics
	_spdiag=tfs->make<TTree>("spdiag","stereo position diagnostics");
	_spdiag->Branch("shphi",&_shphi,"shphi/F");
	_spdiag->Branch("shrho",&_shrho,"shrho/F");
	_spdiag->Branch("stphi",&_stphi,"stphi/F");
	_spdiag->Branch("strho",&_strho,"strho/F");
	_spdiag->Branch("stereo",&_stereo,"stereo/I");
	_spdiag->Branch("chisq",&_schi2,"chisq/F");
	_spdiag->Branch("mvaout",&_smvaout,"mvaout/F");
	_spdiag->Branch("dist",&_sdist,"dist/F");
	_spdiag->Branch("dz",&_sdz,"dz/F");
	_spdiag->Branch("ddot",&_sddot,"ddot/F");
	_spdiag->Branch("mcrel",&_mcr,"mcr/I");
	_spdiag->Branch("mcpdg",&_mcpdg,"mcpdg/I");
	_spdiag->Branch("mcgen",&_mcgen,"mcgen/I");
	_spdiag->Branch("mcproc",&_mcproc,"mcproc/I");
	_spdiag->Branch("mcshphi",&_mcshphi,"mcshphi/F");
	_spdiag->Branch("mcshrho",&_mcshrho,"mcshrho/F");
	_spdiag->Branch("fs",&_sfs,"fs/I");
	if(_diag > 2){
	  _sdiag=tfs->make<TTree>("sdiag","stereo diagnostics");
	  _sdiag->Branch("fs",&_fs,"fs/I");
	  _sdiag->Branch("de",&_de,"de/F");
	  _sdiag->Branch("dt",&_dt,"dt/F");
	  _sdiag->Branch("dz",&_dz,"dz/F");
	  _sdiag->Branch("rho",&_rho,"rho/F");
	  _sdiag->Branch("dist",&_dist,"dist/F");
	  _sdiag->Branch("dperp",&_dperp,"dperp/F");
	  _sdiag->Branch("dl1",&_dl1,"dl1/F");
	  _sdiag->Branch("dl2",&_dl2,"dl2/F");
	  _sdiag->Branch("dc1",&_dc1,"dc1/F");
	  _sdiag->Branch("dc2",&_dc2,"dc2/F");
	  _sdiag->Branch("chi2",&_chi2,"chi2/F");
	  _sdiag->Branch("ndof",&_ndof,"ndof/F");
	  _sdiag->Branch("tdiv1",&_tdiv1,"tdiv1/B");
	  _sdiag->Branch("tdiv2",&_tdiv2,"tdiv2/B");
	  _sdiag->Branch("mvaout",&_mvaout,"mvaout/F");
	  _sdiag->Branch("ddot",&_ddot,"ddot/F");
	  _sdiag->Branch("mcrel",&_mcrel,"mcrel/I");
	  _sdiag->Branch("mcdist",&_mcdist,"mcdist/F");
	  _sdiag->Branch("mcperp",&_mcperp,"mcperp/F");
	}
      }
    }
  }

  void StereoHitDiag::analyze(const art::Event& evt ) {
    const Tracker& tracker = getTrackerOrThrow();
    const TTracker& tt = dynamic_cast<const TTracker&>(tracker);
    findData(evt);
    // draw the stations if requested
    static bool first(false);
    if(_drawstations && first){
      first = false;
      drawStations();
    }
// loop over stereo hits
    for(auto sth : *_stcol) {
      StrawHit const& sh1 = _shcol->at(sth.hitIndex1());
      StrawHit const& sh2 = _shcol->at(sth.hitIndex2());
      StrawHitPosition const& shp1 = _shpcol->at(sth.hitIndex1());
      StrawHitPosition const& shp2 = _shpcol->at(sth.hitIndex2());
      Straw const& straw1 = tt.getStraw(sh1.strawIndex());
      Straw const& straw2 = tt.getStraw(sh2.strawIndex());
    // straw pair diagnostics
      _deltat->Fill(sth.dt());
      double de = fabs(sh1.energyDep()-sh2.energyDep())/sth.energy();
      _deltaE->Fill(de);
      CLHEP::Hep3Vector dpos = shp1.pos() -shp2.pos();
      _deltaz->Fill(dpos.z());
      _sep->Fill(sth.panelSeparation());
      _ddoth->Fill(sth.wdot());
      _dperph->Fill(dpos.perp());
      double dl1 = straw1.getDetail().activeHalfLength()-fabs(sth.wdist1());
      double dl2 = straw2.getDetail().activeHalfLength()-fabs(sth.wdist2());
      _dL->Fill(dl1);
      _dL->Fill(dl2);
      _mva->Fill(sth.mvaout());
      if(_diag > 2){
	_fs = sth.panelSeparation();
	_de = de;
	_dt = sth.dt();
	_dz = dpos.z();
	_rho = sth.pos().perp();
	_dist = dpos.mag();
	_dperp = dpos.perp();
	_dl1 = dl1;
	_dl2 = dl2;
	_dc1 = (shp1.pos()-sth.pos()).perp()/shp1.posRes(StrawHitPosition::wire);
	_dc2 = (shp2.pos()-sth.pos()).perp()/shp2.posRes(StrawHitPosition::wire);
	_tdiv1 = shp1.flag().hasAllProperties(StrawHitFlag::tdiv);
	_tdiv2 = shp2.flag().hasAllProperties(StrawHitFlag::tdiv);
	_chi2 = sth.chisq();
	_ndof = 0; if(_tdiv1)_ndof++; if(_tdiv2) _ndof++;
	_mvaout = sth.mvaout();
	_ddot = sth.wdot();
	_mcdist = _mcperp = -1.0;
	if(_mcdiag){
	  StrawDigiMC const& mcd1 = _mcdigis->at(sth.hitIndex1());
	  StrawDigiMC const& mcd2 = _mcdigis->at(sth.hitIndex2());
	  _mcrel = MCRelationship::relationship(mcd1,mcd2);
	  if(mcd1.stepPointMC(TrkTypes::cal).isNonnull() &&
	      mcd2.stepPointMC(TrkTypes::cal).isNonnull() ){
	    CLHEP::Hep3Vector mcsep = mcd1.stepPointMC(TrkTypes::cal)->position() -
	      mcd2.stepPointMC(TrkTypes::cal)->position();
	    _mcdist = mcsep.mag();
	    _mcperp = mcsep.perp();
	  }
	}
	_sdiag->Fill();
      }
    }
    // loop over straw hit positions
    if(_diag > 1){
      for (size_t ish = 0; ish < _shpcol->size(); ++ish) {
	StrawHitPosition const& sthp = _sthpcol->at(ish);
	StrawHitPosition const& shp = _shpcol->at(ish);
	CLHEP::Hep3Vector shpos = shp.pos(); 
	CLHEP::Hep3Vector sthpos = sthp.pos(); 
	_shphi = shpos.phi();
	_shrho = shpos.perp();
	_stphi = sthpos.phi();
	_strho = sthpos.perp();
	_stereo = sthp.flag().hasAllProperties(StrawHitFlag::stereo);
	_mcr = _sfs = -1;
	_schi2 = _sdist = _sdz = _sddot = -1.0;
	if(sthp.stereoHitIndex() > 0){
	  StereoHit const& sth = _stcol->at(sthp.stereoHitIndex());
	  _schi2 = sth.chisq();
	  _smvaout = sth.mvaout();
	  _sfs = sth.panelSeparation();
	  StrawHitPosition const& shp1 = _shpcol->at(sth.hitIndex1());
	  StrawHitPosition const& shp2 = _shpcol->at(sth.hitIndex2());
	  _sdist = (shp1.pos()-shp2.pos()).mag();
	  _sddot = sth.wdot();
	  _sdz = fabs(shp1.pos().z()-shp2.pos().z());
	  StrawDigiMC const& mcd1 = _mcdigis->at(sth.hitIndex1());
	  StrawDigiMC const& mcd2 = _mcdigis->at(sth.hitIndex2());
	  _mcr = MCRelationship::relationship(mcd1,mcd2);
	}
	StrawDigiMC const& mcd = _mcdigis->at(ish);
	StrawEnd itdc;
	_mcshphi = _mcshrho = -1000.0;
	_mcpdg = _mcproc = _mcgen = 0;
	if(mcd.stepPointMC(itdc).isNonnull() ){
	  _mcshphi = mcd.stepPointMC(itdc)->position().phi();
	  _mcshrho = mcd.stepPointMC(itdc)->position().perp();
	  if(mcd.stepPointMC(itdc)->simParticle().isNonnull()){
	    _mcpdg = mcd.stepPointMC(itdc)->simParticle()->pdgId();
	    _mcproc = mcd.stepPointMC(itdc)->simParticle()->creationCode();
	    if(mcd.stepPointMC(itdc)->simParticle()->genParticle().isNonnull()){
	      _mcgen = mcd.stepPointMC(itdc)->simParticle()->genParticle()->generatorId().id();
	    }
	  }
	}
	_spdiag->Fill();
      }
    }
  }


  bool StereoHitDiag::findData(const art::Event& evt){
    _shcol = 0; _stcol = 0; _sthpcol = 0; _shpcol = 0;  _mcdigis = 0;
    auto shH = evt.getValidHandle<StrawHitCollection>(_shTag);
    _shcol = shH.product();
    auto shpH = evt.getValidHandle<StrawHitPositionCollection>(_shpTag);
    _shpcol = shpH.product();
    auto sthpH = evt.getValidHandle<StrawHitPositionCollection>(_sthpTag);
    _sthpcol = sthpH.product();
    auto stH = evt.getValidHandle<StereoHitCollection>(_stTag);
    _stcol = stH.product();
    if(_mcdiag){
      auto mcdH = evt.getValidHandle<StrawDigiMCCollection>(_mcdigisTag);
      _mcdigis = mcdH.product();
    }
    return _shcol != 0 && _stcol != 0 && _shpcol != 0 && (_mcdigis != 0 || !_mcdiag);
  }

  void StereoHitDiag::drawStations() {
    // station layout diagnostics
    const Tracker& tracker = getTrackerOrThrow();
    const TTracker& tt = dynamic_cast<const TTracker&>(tracker);

    art::ServiceHandle<art::TFileService> tfs;
    unsigned nsta = tt.nPlanes()/2;
    for(unsigned ista=0;ista<nsta;++ista){
      char name[100];
      snprintf(name,100,"station%i",ista);
      _stations.push_back( tfs->make<TH2F>(name,name,100,-700,700,100,-700,700));
      _stations[ista]->SetStats(false);
      TList* flist = _stations[ista]->GetListOfFunctions();
      TLegend* sleg = new TLegend(0.1,0.6,0.3,0.9);
      flist->Add(sleg);
      for(int iplane=0;iplane<2;++iplane){
	const Plane& pln = tt.getPlane(2*ista+iplane);
	const vector<Panel>& panels = pln.getPanels();
	for(size_t ipnl=0;ipnl<panels.size();++ipnl){
	  int iface = ipnl%2;
	  const Panel& pnl = panels[ipnl];
	  CLHEP::Hep3Vector spos = pnl.straw0MidPoint();
	  CLHEP::Hep3Vector sdir = pnl.straw0Direction();
	  CLHEP::Hep3Vector end0 = spos - 100.0*sdir;
	  CLHEP::Hep3Vector end1 = spos + 100.0*sdir;
	  TLine* sline = new TLine(end0.x(),end0.y(),end1.x(),end1.y());
	  sline->SetLineColor(ipnl+1);
	  sline->SetLineStyle(2*iplane+iface+1);
	  flist->Add(sline);
	  TMarker* smark = new TMarker(end0.x(),end0.y(),8);
	  smark->SetMarkerColor(ipnl+1);
	  smark->SetMarkerSize(2);
	  flist->Add(smark);
	  char label[80];
	  snprintf(label,80,"pln %i pnl %i",iplane,(int)ipnl);
	  sleg->AddEntry(sline,label,"l");
	}
      }
    }
  }
}

// Part of the magic that makes this class a module.
using mu2e::StereoHitDiag;
DEFINE_ART_MODULE(StereoHitDiag);

