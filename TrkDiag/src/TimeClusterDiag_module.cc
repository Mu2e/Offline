//
// TTracker time cluster diagnostics
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
#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
// tracking
#include "TrkReco/inc/TrkTimeCalculator.hh"
#include "TrkReco/inc/TrkUtilities.hh"
// TrkDiag
#include "TrkDiag/inc/TimeClusterInfo.hh"
#include "TrkDiag/inc/TrkMCTools.hh"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
// root
#include "TH1F.h"
#include "TTree.h"
#include "TMarker.h"
// C++
#include <functional>
#include <iostream>
#include <algorithm>
using namespace std; 
using CLHEP::Hep3Vector;

namespace mu2e {
// ttree structs
  class TimeClusterDiag : public art::EDAnalyzer {
  public:

    explicit TimeClusterDiag(fhicl::ParameterSet const& pset);
    virtual ~TimeClusterDiag();

    virtual void beginJob();

    // This is called for each event.
    virtual void analyze(art::Event const& e);

  private:

    int           _diag;
    bool	  _mcdiag;
    bool	  _plotts;
    unsigned _minnce; // minimum # CE hits to make plots
    double _ccmine; // min calo cluster energy to plot
    double _ccwt;

    // event object Tags
    art::InputTag   _shTag;
    art::InputTag   _shpTag;
    art::InputTag   _shfTag;
    art::InputTag   _tcTag;
    art::InputTag   _ccTag;
    art::InputTag   _mcdigisTag;
    // hit selectors
    StrawHitFlag  _hsel, _tcsel, _hbkg;

    // mva stuff
    MVATools                              _peakMVA; // MVA for peak cleaning
    // time spectrum parameters
    double        _tmin;
    double        _tmax;
    double        _tbin;
    unsigned      _nbins;
    // cache of event objects
    const StrawHitCollection*             _shcol;
    const StrawHitFlagCollection*         _shfcol;
    const StrawHitPositionCollection*     _shpcol;
    const TimeClusterCollection*	  _tccol;
    const CaloClusterCollection*	  _cccol;
    const StrawDigiMCCollection*          _mcdigis;
// TTree variables
    TTree*                                _tcdiag;
    // event number
    Int_t      _iev;
    Int_t     _ntc; // # clusters/event
    TimeClusterInfo		  _besttc;  // info about best time cluster (most CE hits)
    MCClusterInfo		  _ceclust; // info about 'cluster' of MC CE hits
    vector<TimeClusterHitInfo>	  _tchinfo; // info about hits of best time cluster
    vector<TimeClusterInfo>	  _alltc; // info about all TimeClusters in the event, sorted by # CE
    // time offsets
    SimParticleTimeOffset _toff;
    //t0 calculator
    TrkTimeCalculator _ttcalc;
// helper functions
    void createDiagnostics();
    bool findData(art::Event const& evt);
    void plotTimeSpectra ();
    void fillClusterInfo ();
    void fillClusterInfo (TimeCluster const& tc,TimeClusterInfo& tcinfo);
    void fillClusterHitInfo (TimeCluster const& besttc); 
    void fillCECluster();
  };

  TimeClusterDiag::~TimeClusterDiag() {
  }
  
  TimeClusterDiag::TimeClusterDiag(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _diag		(pset.get<int>("diagLevel",1)),
    _mcdiag		(pset.get<bool>("MCdiag",true)),
    _plotts		(pset.get<bool>("PlotTimeSpectra",false)),
    _minnce		(pset.get<unsigned>("MinimumCEHits",15)),
    _ccmine		(pset.get<double>("CaloClusteriMinE",50.0)),
    _ccwt		(pset.get<double>("CaloClusterWeight",10.0)),
    _shTag		(pset.get<art::InputTag>("StrawHitCollection","makeSH")),
    _shpTag		(pset.get<art::InputTag>("StrawHitPositionCollection","MakeStereoHits")),
    _shfTag		(pset.get<art::InputTag>("StrawHitFlagCollection","TimeClusterFinder")),
    _tcTag		(pset.get<art::InputTag>("TimeClusterCollection","TimeClusterFinder")),
    _ccTag              (pset.get<art::InputTag>("caloClusterModuleLabel","MakeCaloCluster")),
    _mcdigisTag		(pset.get<art::InputTag>("StrawDigiMCCollection","makeSH")),
    _hsel		(pset.get<std::vector<std::string> >("HitSelectionBits",vector<string>{"EnergySelection","TimeSelection","RadiusSelection"})),
    _hbkg		(pset.get<vector<string> >("HitBackgroundBits",vector<string>{"Background"})),
    _peakMVA		(pset.get<fhicl::ParameterSet>("PeakCleanMVA",fhicl::ParameterSet())),
    _tmin		(pset.get<double>("tmin",500.0)),
    _tmax		(pset.get<double>("tmax",1700.0)),
   _tbin		(pset.get<double>("tbin",15.0)),
    _toff(pset.get<fhicl::ParameterSet>("TimeOffsets")),
    _ttcalc            (pset.get<fhicl::ParameterSet>("T0Calculator",fhicl::ParameterSet()))
    {
    // set # bins for time spectrum plot
    _nbins = (unsigned)rint((_tmax-_tmin)/_tbin);
    _tcsel = StrawHitFlag(StrawHitFlag::tclust); 
  }

  void TimeClusterDiag::beginJob(){
  // initialize MVA: this is used just for diagnostics
    _peakMVA.initMVA();
    createDiagnostics();
  }

  void TimeClusterDiag::analyze(art::Event const& event ) {
    _iev=event.id().event();
    // find the data
    if(!findData(event) || _shcol->size() != _shfcol->size() 
    || _shcol->size() != _shpcol->size()
    || (_mcdigis !=0 && _shcol->size() != _mcdigis->size())){
      throw cet::exception("RECO")<<"mu2e::TimeClusterDiag: data missing or inconsistent"<< endl;
      return;
    }
    _ntc = _tccol->size();
    // fill MC info
    if(_mcdiag){
      fillCECluster();
    }
    // fill info for all TimeClusters
    fillClusterInfo();
    if(_alltc.size() > 0)
      // best is defined as having the most CE hits
      _besttc = _alltc[0];
    else
      _besttc.reset();
    if(_diag > 1 && _besttc._tcindex >= 0 && _besttc._tcindex < static_cast<int>(_tccol->size())){
      _tchinfo.clear();
      fillClusterHitInfo(_tccol->at(_besttc._tcindex));
    }
    // fill the tree
    _tcdiag->Fill();
    // if requested, plot time spectra for this event
    if(_plotts){
      unsigned nce(0);
      for(auto const& mcdigi : *_mcdigis) {
	if(TrkMCTools::CEDigi(mcdigi))++nce;
      }
      if (nce >= _minnce) 
	plotTimeSpectra();
    }
  } // end analyze

  // find the input data objects
  bool TimeClusterDiag::findData(art::Event const& evt){
    _shcol = 0; _shfcol = 0; _shpcol = 0; _tccol = 0; _mcdigis = 0;
    auto shH = evt.getValidHandle<StrawHitCollection>(_shTag);
    _shcol = shH.product();
    auto shpH = evt.getValidHandle<StrawHitPositionCollection>(_shpTag);
    _shpcol = shpH.product();
    auto shfH = evt.getValidHandle<StrawHitFlagCollection>(_shfTag);
    _shfcol = shfH.product();
    auto tcH = evt.getValidHandle<TimeClusterCollection>(_tcTag);
    _tccol = tcH.product();
    if(_mcdiag){
      auto mcdH = evt.getValidHandle<StrawDigiMCCollection>(_mcdigisTag);
      _mcdigis = mcdH.product();
      _toff.updateMap(evt);
    }
    // calorimeter data may or may not be present
    art::Handle<CaloClusterCollection> ccH;
    evt.getByLabel<CaloClusterCollection>(_ccTag,ccH);
    if(ccH.isValid())
      _cccol = ccH.product();

    return _shcol != 0 && _shfcol != 0 && _shpcol != 0 && _tccol != 0 && (_mcdigis != 0 || !_mcdiag);
  }

  void TimeClusterDiag::fillCECluster() {
    _ceclust.reset();
  // loop over all straw hits 
    unsigned nstrs = _shcol->size();
    Hep3Vector cpos;
    for(unsigned istr=0; istr<nstrs;++istr){
      StrawDigiMC const& mcdigi = _mcdigis->at(istr);
      StrawHitFlag const& shflag =_shfcol->at(istr);
      if(TrkMCTools::CEDigi(mcdigi)){
	++_ceclust._nce;
	StrawDigi::TDCChannel itdc = StrawDigi::zero;
	if(!mcdigi.hasTDC(itdc))itdc = StrawDigi::one;
	cpos += mcdigi.clusterPosition(itdc).vect();
	_ceclust._time += mcdigi.clusterPosition(itdc).t();
	bool selected = shflag.hasAllProperties(_hsel) && !shflag.hasAnyProperty(_hbkg);
	bool tclust = shflag.hasAllProperties(_tcsel);
	if(selected)++_ceclust._ncesel;
	if(tclust)++_ceclust._nceclust;
      }
    }
    if(_ceclust._nce > 0){
      _ceclust._pos = cpos/_ceclust._nce;
      _ceclust._time /=_ceclust._nce;
    }
// 2nd pass to get extents {
    double cphi = cpos.phi();
    for(unsigned istr=0; istr<nstrs;++istr){
      StrawDigiMC const& mcdigi = _mcdigis->at(istr);
      if(TrkMCTools::CEDigi(mcdigi)){
	StrawDigi::TDCChannel itdc = StrawDigi::zero;
	if(!mcdigi.hasTDC(itdc))itdc = StrawDigi::one;
	Hep3Vector hpos = mcdigi.clusterPosition(itdc).vect();
	double hphi = hpos.phi();
	float hrho = hpos.perp();
	float dphi = fabs(Angles::deltaPhi(hphi,cphi));
	_ceclust._maxdphi = std::max(dphi,_ceclust._maxdphi);
	_ceclust._minrho = std::min(hrho,_ceclust._minrho);
	_ceclust._maxrho = std::max(hrho,_ceclust._maxrho);
      }
    }
  }

  void TimeClusterDiag::fillClusterHitInfo(TimeCluster const& tc) {
    for (auto const& idx : tc._strawHitIdxs) {
      TimeClusterHitInfo tchi;
      size_t ish = idx;
      tchi._time = _shcol->at(ish).time();
      tchi._dt = _ttcalc.strawHitTime(_shcol->at(ish),_shpcol->at(ish)) -tc._t0._t0;
      Hep3Vector const& pos = _shpcol->at(ish).pos();
      double phi = pos.phi();
      tchi._dphi = Angles::deltaPhi(phi,tc._pos.phi());
      tchi._rho = pos.perp();
      tchi._z = pos.z();
// compute MVA
      std::vector<Double_t> pars(3);
      pars[0] = tchi._dt;
      pars[1] = tchi._dphi;
      pars[2] = tchi._rho;
      tchi._mva = _peakMVA.evalMVA(pars);
// MC truth
      if(_mcdiag){
	StrawDigiMC const& mcdigi = _mcdigis->at(idx);
	StrawDigi::TDCChannel itdc = StrawDigi::zero;
	if(!mcdigi.hasTDC(itdc))itdc= StrawDigi::one;
	tchi._mctime = _toff.timeWithOffsetsApplied( *mcdigi.stepPointMC(itdc));
	tchi._mcmom = mcdigi.stepPointMC(itdc)->momentum().mag();
	art::Ptr<SimParticle> sp;
	if(TrkMCTools::simParticle(sp,mcdigi) > 0){
	  tchi._mcpdg = sp->pdgId();
	  tchi._mcproc = sp->creationCode();
	  if(sp->genParticle().isNonnull())
	    tchi._mcgen = sp->genParticle()->generatorId().id();
	}
      }
      _tchinfo.push_back(tchi);
    }
  }

  void TimeClusterDiag::plotTimeSpectra() {
    art::ServiceHandle<art::TFileService> tfs;
    TH1F *rtsp, *stsp, *ctsp, *actsp, *sctsp, *cctsp, *catsp;

    char rsname[100];
    char ssname[100];
    char csname[100];
    char acsname[100];
    char scsname[100];
    char ccsname[100];
    char casname[100];

    snprintf(rsname,100,"rawtspectrum%i",_iev);
    snprintf(ssname,100,"seltspectrum%i",_iev);
    snprintf(csname,100,"clusttspectrum%i",_iev);
    snprintf(acsname,100,"allconvtspectrum%i",_iev);
    snprintf(scsname,100,"selconvtspectrum%i",_iev);
    snprintf(ccsname,100,"clustconvtspectrum%i",_iev);
    snprintf(casname,100,"calotspectrum%i",_iev);
    
    rtsp = tfs->make<TH1F>(rsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    rtsp->SetLineColor(kCyan);
    stsp = tfs->make<TH1F>(ssname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    stsp->SetLineColor(kGreen);
    ctsp = tfs->make<TH1F>(csname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    ctsp->SetLineColor(kBlue);
    actsp = tfs->make<TH1F>(acsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    actsp->SetLineColor(kRed);
    sctsp = tfs->make<TH1F>(scsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    sctsp->SetLineColor(kYellow);
    sctsp->SetFillColor(kYellow);
    sctsp->SetFillStyle(3001);
    cctsp = tfs->make<TH1F>(ccsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    cctsp->SetLineColor(kRed);
    catsp = tfs->make<TH1F>(casname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    catsp->SetLineColor(kMagenta);

    unsigned nstrs = _shcol->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      double time = _ttcalc.strawHitTime(_shcol->at(istr),_shpcol->at(istr));
      bool conversion(false);
      if(_mcdigis != 0) {
	StrawDigiMC const& mcdigi = _mcdigis->at(istr);
	conversion = TrkMCTools::CEDigi(mcdigi);
      }
      bool selected = _shfcol->at(istr).hasAllProperties(_hsel) && !_shfcol->at(istr).hasAnyProperty(_hbkg);
      bool tclust = _shfcol->at(istr).hasAllProperties(_tcsel);
      // fill plots
      rtsp->Fill(time);
      if(selected)stsp->Fill(time);
      if(tclust)ctsp->Fill(time);
      if(conversion){
	actsp->Fill(time);
	if(selected)sctsp->Fill(time);
	if(tclust)cctsp->Fill(time);
      }
    }
    // clusters if available
    if(_cccol != 0){
      for(auto cc : *_cccol) {
	if(cc.energyDep() > _ccmine)
	  catsp->Fill( _ttcalc.caloClusterTime(cc),_ccwt);
      }
    }
    // plot time cluster times
    TList* flist = ctsp->GetListOfFunctions();
    for(auto itclust=_tccol->begin(); itclust!=_tccol->end(); ++itclust){
    // change marker if calo cluster is present
      int imark = 23;
      if(itclust->_caloCluster.isNonnull()) imark = 34;
      TMarker* smark = new TMarker(itclust->_t0._t0,ctsp->GetBinContent(ctsp->FindBin(itclust->_t0._t0)),imark);
      smark->SetMarkerColor(kRed);
      smark->SetMarkerSize(1.5);
      flist->Add(smark);
    }
  }

  void TimeClusterDiag::fillClusterInfo () {
    _alltc.clear();
    for(size_t ic = 0;ic< _tccol->size(); ++ ic){
      TimeCluster const& tc = _tccol->at(ic);
      TimeClusterInfo tcinfo;
      tcinfo._tcindex = ic;
      fillClusterInfo(tc,tcinfo);
      _alltc.push_back(tcinfo);
    }
// sort by the # of CE
    if(_mcdiag){
      static NCEComp comp;
      std::sort(_alltc.begin(),_alltc.end(),comp);
    }
  }
    
  void TimeClusterDiag::fillClusterInfo (TimeCluster const& tp,TimeClusterInfo& tcinfo) {
    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    const Calorimeter* calo = ch.get();
// simple entries
    tcinfo._nhits = tp._strawHitIdxs.size();
    tcinfo._time  = tp._t0._t0;
    tcinfo._terr  = tp._t0._t0err;
    tcinfo._pos	  = tp._pos;
    // calo info if available
    if(tp._caloCluster.isNonnull()){
      tcinfo._ecalo = tp._caloCluster->energyDep();
      tcinfo._tcalo = tp._caloCluster->time();
      tcinfo._dtcalo = _ttcalc.caloClusterTime(*tp._caloCluster) - tp._t0._t0;
      // calculate the cluster position.  Currently the Z is in disk coordinates and must be translated, FIXME!
      Hep3Vector cog = calo->geomUtil().mu2eToTracker(calo->geomUtil().diskFFToMu2e(tp._caloCluster->diskId(),tp._caloCluster->cog3Vector()));
      tcinfo._cog = cog;
    }
    // hit summary
    tcinfo._minhtime = 1.0e9;
    tcinfo._maxhtime = 0.0; 
    tcinfo._ncehits = 0;
    // count CE hits 
    for (auto const& idx : tp._strawHitIdxs) {
      float ht = _ttcalc.strawHitTime(_shcol->at(idx),_shpcol->at(idx));
      tcinfo._maxhtime = std::max( ht , tcinfo._maxhtime);
      tcinfo._minhtime = std::min( ht , tcinfo._minhtime);
      if(_mcdiag){
      // mc truth info
	StrawDigiMC const& mcdigi = _mcdigis->at(idx);
	bool conversion = TrkMCTools::CEDigi(mcdigi);
	if(conversion)++tcinfo._ncehits;
      }
    }
    // look for overlaps
    for(size_t ic = 0;ic< _tccol->size(); ++ ic){
      if((int)ic != tcinfo._tcindex) // don't count myself
	tcinfo._maxover = max(tcinfo._maxover,(Float_t)TrkUtilities::overlap(tp,_tccol->at(ic)));
    }
  }

  void TimeClusterDiag::createDiagnostics() {
    art::ServiceHandle<art::TFileService> tfs;
   // time peak diagnostics
    _tcdiag=tfs->make<TTree>("tcdiag","time cluster diagnostics");
    // event info should use the TrkAna struct, FIXME!!
    _tcdiag->Branch("iev",&_iev,"iev/I");
    _tcdiag->Branch("ntc",&_ntc,"ntc/I");
    _tcdiag->Branch("besttc",&_besttc,TimeClusterInfo::leafnames().c_str());
    if(_mcdiag){
      _tcdiag->Branch("ceclust",&_ceclust,MCClusterInfo::leafnames().c_str());
    }
    if(_diag > 1) _tcdiag->Branch("tchinfo",&_tchinfo);
    if(_diag > 2) _tcdiag->Branch("alltc",&_alltc);
  }




}  // end namespace mu2e

using mu2e::TimeClusterDiag;
DEFINE_ART_MODULE(TimeClusterDiag);
