//
// Tracker time cluster diagnostics
//
// Original author D. Brown and G. Tassielli
//

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art_root_io/TFileService.h"
// mu2e
#include "Offline/GeneralUtilities/inc/Angles.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
#include "Offline/MCDataProducts/inc/MCRelationship.hh"
// tracking
#include "Offline/TrkReco/inc/TrkTimeCalculator.hh"
#include "Offline/TrkReco/inc/TrkUtilities.hh"
// TrkDiag
#include "Offline/TrkDiag/inc/TimeClusterInfo.hh"
#include "Offline/TrkDiag/inc/TrkMCTools.hh"
// data
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
// root
#include "TH1F.h"
#include "TTree.h"
#include "TMarker.h"
// C++
#include <functional>
#include <iostream>
#include <algorithm>
using namespace std;

namespace mu2e {
  // ttree structs
  class TimeClusterDiag : public art::EDAnalyzer {
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<int> diag{ Name("diagLevel"), Comment("Diag"),0 };
      fhicl::Atom<bool> mcdiag{ Name("MCDiag"), Comment("MonteCarlo Diag"), true};
      fhicl::Atom<bool> plotTS{ Name("PlotTimeSpectra"), Comment("Plot individual event time spectra"), false };
      fhicl::Atom<int> minNCE{ Name("MinimumNCEHits"), Comment("Mininum # Ce hits to plot event time spectra"), 15 };

      fhicl::Sequence<std::string> hsel { Name("HitSelectionBits"), Comment("Good hit bits") };
      fhicl::Sequence<std::string> hbkg { Name("HitBackgroundBits"), Comment("Background hit bits") };
      fhicl::Atom<double> tmin { Name("tmin"), Comment("Histogram lower time"),500.0 };
      fhicl::Atom<double> tmax { Name("tmax"), Comment("Histogram upper time"),1700.0};
      fhicl::Atom<double> tbin { Name("tbin"), Comment("Histogram time bin"),15.0};
      fhicl::Atom<double> pitch { Name("AveragePitch"), Comment("Average Pitch for time correctin"),0.6};

      fhicl::Atom<art::InputTag> ComboHitCollection{   Name("ComboHitCollection"),   Comment("ComboHit collection name") };
      fhicl::Atom<art::InputTag> StrawDigiMCCollection{   Name("StrawDigiMCCollection"),   Comment("StrawDigiMC collection name") };
      fhicl::Atom<art::InputTag> TimeClusterCollection{   Name("TimeClusterCollection"),   Comment("TimeCluster collection name") };
      fhicl::Atom<art::InputTag> MCPrimary{ Name("MCPrimary"),Comment("MC Primary Particle") };
      fhicl::Atom<art::InputTag> VDMC{ Name("VDStepPointMCCollection"),Comment("Virtual Detector StepPointMC collection name") };
      fhicl::Table<TrkTimeCalculator::Config> ttcalc {Name("T0Calculator"),           Comment("TimeTracker calculator config") };
    };

    // MC matching
    struct spcount {
      spcount() : _count(0) {}
      spcount(art::Ptr<SimParticle> const& spp) : _spp(spp), _count(1) {}
      void append(art::Ptr<SimParticle> const& sp) { if(sp == _spp)++_count; }
      bool operator ==(art::Ptr<SimParticle> const& sp) const { return _spp == sp; }
      art::Ptr<SimParticle> _spp;
      unsigned _count;
    };
    struct spcountcomp {
      bool operator() (spcount a, spcount b) { return a._count > b._count; }
    };

    public:
    explicit TimeClusterDiag(const art::EDAnalyzer::Table<Config>& config);
    virtual ~TimeClusterDiag();
    virtual void beginJob();
    virtual void analyze(art::Event const& e);

    private:

    int  _diag;
    bool _mcdiag;
    bool    _plotts;
    // event objects
    art::ProductToken<ComboHitCollection> _chToken;
    art::ProductToken<TimeClusterCollection> _tcToken;
    art::ProductToken<StrawDigiMCCollection> _mcdigisToken;
    art::ProductToken<PrimaryParticle> _mcprimaryToken;
    art::ProductToken<StepPointMCCollection> _vdmcstepsToken;
    // hit selectors
    StrawHitFlag  _hsel, _tcsel, _hbkg;

    // time spectrum parameters
    double        _tmin;
    double        _tmax;
    double        _tbin;
    unsigned      _nbins;
    // cache of event objects
    const ComboHitCollection*    _chcol;
    const TimeClusterCollection* _tccol;
    const StrawDigiMCCollection* _mcdigis;
    const StepPointMCCollection* _vdmcsteps;
    const PrimaryParticle *_mcprimary;
    StrawHitFlag _cesel; // flag bit for Ce (from truth)
    // TTree variables
    TTree*  _tcdiag;
    // event number
    Int_t     _iev;
    Int_t     _ntc; // # clusters/event
    Float_t   _mcmidt0;
    TimeClusterInfo     _besttc;  // info about best time cluster (most CE hits)
    MCClusterInfo     _ceclust; // info about 'cluster' of MC CE hits
    vector<TimeClusterHitInfo>    _tchinfo; // info about hits of best time cluster
    vector<TimeClusterInfo>   _alltc; // info about all TimeClusters in the event, sorted by # CE
    float   _pitch; // average helix pitch (= dz/dflight, =sin(lambda))
    StrawHitFlagCollection _shfcol;
    //t0 calculator
    TrkTimeCalculator _ttcalc;
    unsigned _minnce;
    // helper functions
    void createDiagnostics();
    bool findData(art::Event const& evt);
    void setFlags(art::Event const& evt);
    void plotTimeSpectra ();
    void fillClusterInfo (std::vector<spcount> const& primaries);
    void fillClusterInfo (TimeCluster const& tc,spcount const& primary, TimeClusterInfo& tcinfo);
    void fillClusterHitInfo (TimeCluster const& besttc,
        art::Ptr<SimParticle> const& primary, art::Event const& evt);
    void fillCECluster();
    void findPrimaries (art::Event const& evt, std::vector<spcount>& primaries);
    void findPrimary(art::Event const& evt, TimeCluster const& tc, spcount& primary );
  };

  TimeClusterDiag::~TimeClusterDiag() {
  }

  TimeClusterDiag::TimeClusterDiag(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer(config),
    _diag   (config().diag()),
    _mcdiag   (config().mcdiag()),
    _plotts   (config().plotTS()),
    _chToken{ consumes<ComboHitCollection>(config().ComboHitCollection() ) },
    _tcToken{ consumes<TimeClusterCollection>(config().TimeClusterCollection() ) },
    _mcdigisToken{ consumes<StrawDigiMCCollection>(config().StrawDigiMCCollection() ) },
    _mcprimaryToken{ consumes<PrimaryParticle>(config().MCPrimary() ) },
    _vdmcstepsToken{ consumes<StepPointMCCollection>(config().VDMC() ) },
    _hsel   (config().hsel()),
    _hbkg   (config().hbkg()),
    _tmin   (config().tmin()),
    _tmax   (config().tmax()),
    _tbin    (config().tbin()),
    _pitch   (config().pitch()),
    _ttcalc  (config().ttcalc()),
    _minnce (config().minNCE()) {
      // set # bins for time spectrum plot
      _nbins = (unsigned)rint((_tmax-_tmin)/_tbin);
      _tcsel = StrawHitFlag(StrawHitFlag::tclust);
      _cesel = StrawHitFlag(StrawHitFlag::track);
    }

  void TimeClusterDiag::beginJob(){
    createDiagnostics();
  }

  void TimeClusterDiag::analyze(art::Event const& event ) {
    _iev=event.id().event();
    // find the data
    if(!findData(event) )
      throw cet::exception("RECO")<<"mu2e::TimeClusterDiag: data missing or inconsistent"<< endl;
    _ntc = _tccol->size();
    //
    setFlags(event);
    // fill MC info
    std::vector<spcount> primaries;
    if(_mcdiag){
      fillCECluster();
      findPrimaries(event, primaries );
      // get _mcmidt0
      art::Ptr<SimParticle> spp;
      _mcmidt0 = -1000;
      // loop over the digis and find the ones that match
      cet::map_vector_key trkid = _mcprimary->primarySimParticles().front()->id();
      auto jmc = _vdmcsteps->end();
      for(auto imc = _vdmcsteps->begin();imc != _vdmcsteps->end(); ++imc ) {
        // find matching steps
        if(  imc->trackId() == trkid &&
            (imc->volumeId() == VirtualDetectorId::TT_Mid || imc->volumeId() == VirtualDetectorId::TT_MidInner)) {
            if(jmc == _vdmcsteps->end() || imc->time() < jmc->time()) // earliest time
              jmc = imc;
        }
      }
      if (jmc != _vdmcsteps->end())
        _mcmidt0 = jmc->time();
    }
    // fill info for all TimeClusters
    fillClusterInfo(primaries);
    if(_alltc.size() > 0)
      // best is defined as having the most CE hits
      _besttc = _alltc[0];
    else
      _besttc.reset();


    if(_diag > 1 && _besttc._tcindex >= 0 && _besttc._tcindex < static_cast<int>(_tccol->size())){
      _tchinfo.clear();
      fillClusterHitInfo(_tccol->at(_besttc._tcindex),
          primaries.at(_besttc._tcindex)._spp,event);
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
     _chcol = 0; _tccol = 0; _mcdigis = 0;
    _vdmcsteps = 0;
    auto chH = evt.getValidHandle(_chToken);
    _chcol = chH.product();
    auto tcH = evt.getValidHandle(_tcToken);
    _tccol = tcH.product();
    if(_mcdiag){
      auto mcdH = evt.getValidHandle(_mcdigisToken);
      _mcdigis = mcdH.product();
      auto mcstepsH = evt.getValidHandle(_vdmcstepsToken);
      _vdmcsteps = mcstepsH.product();
      auto mcpH = evt.getValidHandle(_mcprimaryToken);
      _mcprimary = mcpH.product();
   }
    return _chcol != 0 && _tccol != 0 && ( (_mcdigis != 0 && _vdmcsteps != 0 && _mcprimary != 0) || !_mcdiag);
  }

  void TimeClusterDiag::setFlags(art::Event const& event ) {
    _shfcol.clear();
    _shfcol.reserve(_chcol->size());
    // loop over all hits
    unsigned nch = _chcol->size();
    XYZVectorF cpos;
    for(unsigned ich=0; ich<nch;++ich){
      ComboHit const& ch = _chcol->at(ich);
      _shfcol.push_back(ch.flag());
      if(_mcdiag){
        std::vector<StrawDigiIndex> shids;
        _chcol->fillStrawDigiIndices(ich,shids);
        unsigned nce(0);
        for(auto idigi : shids) {
          StrawDigiMC const& mcdigi = _mcdigis->at(idigi);
          if(TrkMCTools::CEDigi(mcdigi))
            ++nce;
        }
        if(nce/float(shids.size()) >= 0.5)
          _shfcol.back().merge(_cesel);
      }
    }
    // flag hits in clusters
    //
    for(auto const& tc : *_tccol) {
      for (auto ich : tc._strawHitIdxs) {
        _shfcol[ich].merge(_tcsel);
      }
    }
  }

  void TimeClusterDiag::fillCECluster() {
    _ceclust.reset();
    // loop over all straw hits
    unsigned nstrs = _chcol->size();
    unsigned nce(0);
    for(unsigned ich=0; ich<nstrs;++ich){
      ComboHit const& ch = _chcol->at(ich);
      if(_shfcol[ich].hasAllProperties(_cesel)){
        ++nce;
        _ceclust._pos += ch.pos();
        //_ceclust._time +=  _ttcalc.comboHitTime(ch);
        _ceclust._time +=  _ttcalc.comboHitTime(ch,_pitch);
      }
    }

    //_ceclust._pos /= nce;
    //_ceclust._time /= nce;

    if(nce > 0){
      _ceclust._pos /= nce;
      _ceclust._time /= nce;
    }

    // 2nd pass to get extents {
    double cphi = _ceclust._pos.phi();
    double mindphi(1e10), maxdphi(-1e10);
    for(unsigned ich=0; ich<nstrs;++ich){
      ComboHit const& ch = _chcol->at(ich);
      if(_shfcol[ich].hasAllProperties(_cesel)){
        _ceclust._nce += ch.nStrawHits();
        XYZVectorF cpos = ch.pos();

        float hrho = sqrt(cpos.Perp2());
        double hphi = cpos.phi();
        double dphi = Angles::deltaPhi(hphi,cphi);
        maxdphi = std::max(dphi,maxdphi);
        mindphi = std::min(dphi,mindphi);
        _ceclust._minrho = std::min(hrho,_ceclust._minrho);
        _ceclust._maxrho = std::max(hrho,_ceclust._maxrho);
        if(_shfcol[ich].hasAllProperties(_hsel) && !_shfcol[ich].hasAnyProperty(_hbkg))
          _ceclust._ncesel += ch.nStrawHits();
        if(_shfcol[ich].hasAllProperties(_tcsel))
          _ceclust._nceclust += ch.nStrawHits();
      }
    }
    _ceclust._maxdphi = maxdphi-mindphi;
  }

  void TimeClusterDiag::fillClusterHitInfo(TimeCluster const& tc,
      art::Ptr<SimParticle> const& primary,  art::Event const& event) {
    for (auto ich : tc._strawHitIdxs) {
      ComboHit const& ch = _chcol->at(ich);
      TimeClusterHitInfo tchi;
      tchi._time = ch.time();
      tchi._dt = _ttcalc.comboHitTime(ch,_pitch) -tc._t0._t0;
      tchi._wdist = ch.wireDist();
      tchi._werr = ch.wireRes();
      XYZVectorF const& pos = ch.pos();
      double phi = ch.phi();
      tchi._dphi = Angles::deltaPhi(phi,tc._pos.phi());
      tchi._rho = sqrt(pos.Perp2());
      tchi._z = pos.z();
      tchi._edep = ch.energyDep();
      tchi._nsh = ch.nStrawHits();
      tchi._plane = ch.strawId().plane();
      // MC truth
      if(_mcdiag){
        std::vector<StrawDigiIndex> shids;
        _chcol->fillStrawDigiIndices(ich,shids);
        StrawDigiMC const& mcdigi = _mcdigis->at(shids[0]);// FIXME!
        StrawEnd itdc;
        tchi._mctime = mcdigi.strawGasStep(itdc)->time();
        tchi._mcmom = sqrt(mcdigi.strawGasStep(itdc)->momentum().mag2());
        art::Ptr<SimParticle> sp = mcdigi.earlyStrawGasStep()->simParticle();
        tchi._mcpdg = sp->pdgId();
        tchi._mcproc = sp->creationCode();
        if(sp->genParticle().isNonnull())
          tchi._mcgen = sp->genParticle()->generatorId().id();
        if(primary.isNonnull()){
          MCRelationship rel(primary,sp);
          tchi._mcrel = rel.relationship();
        }
      }
      _tchinfo.push_back(tchi);
    }
  }

  void TimeClusterDiag::plotTimeSpectra() {
    art::ServiceHandle<art::TFileService> tfs;
    TH1F *rtsp, *stsp, *ctsp, *actsp, *sctsp, *cctsp;

    char rsname[100];
    char ssname[100];
    char csname[100];
    char acsname[100];
    char scsname[100];
    char ccsname[100];

    snprintf(rsname,100,"rawtspectrum%i",_iev);
    snprintf(ssname,100,"seltspectrum%i",_iev);
    snprintf(csname,100,"clusttspectrum%i",_iev);
    snprintf(acsname,100,"allconvtspectrum%i",_iev);
    snprintf(scsname,100,"selconvtspectrum%i",_iev);
    snprintf(ccsname,100,"clustconvtspectrum%i",_iev);

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

    unsigned nstrs = _chcol->size();
    for(unsigned ich=0; ich<nstrs;++ich){
      ComboHit const& ch = _chcol->at(ich);
      float time = _ttcalc.comboHitTime(ch,_pitch);
      bool conversion = _shfcol[ich].hasAllProperties(_cesel);
      bool selected = _shfcol[ich].hasAllProperties(_hsel) && !_shfcol[ich].hasAnyProperty(_hbkg);
      bool tclust = _shfcol[ich].hasAllProperties(_tcsel);
      // fill plots
      rtsp->Fill(time,ch.nStrawHits());
      if(selected)stsp->Fill(time,ch.nStrawHits());
      if(tclust)ctsp->Fill(time,ch.nStrawHits());
      if(conversion){
        actsp->Fill(time,ch.nStrawHits());
        if(selected)sctsp->Fill(time,ch.nStrawHits());
        if(tclust)cctsp->Fill(time,ch.nStrawHits());
      }
    }
    // plot time cluster times
    TList* flist = ctsp->GetListOfFunctions();
    for(auto itclust=_tccol->begin(); itclust!=_tccol->end(); ++itclust){
      // change marker if calo cluster is present
      int imark = 23;
      if(itclust->_caloCluster.isNonnull()) imark = 34;
      TMarker* smark = new TMarker(itclust->_t0._t0,ctsp->GetMaximum()*1.05,imark);
      smark->SetMarkerColor(kBlue);
      smark->SetMarkerSize(1.5);
      flist->Add(smark);
    }
    // MC true marker
    TList* clist = cctsp->GetListOfFunctions();
    // change marker if calo cluster is present
    TMarker* cmark = new TMarker(_ceclust._time,actsp->GetMaximum()*1.05,29);
    cmark->SetMarkerColor(kRed);
    cmark->SetMarkerSize(1.5);
    clist->Add(cmark);
  }

  void TimeClusterDiag::fillClusterInfo (std::vector<spcount> const& primaries) {
    _alltc.clear();
    for(size_t ic = 0;ic< _tccol->size(); ++ ic){
      TimeCluster const& tc = _tccol->at(ic);
      TimeClusterInfo tcinfo;
      tcinfo._tcindex = ic;
      fillClusterInfo(tc,primaries[ic],tcinfo);
      _alltc.push_back(tcinfo);
    }
    // sort by the # of CE
    if(_mcdiag){
      static NCEComp comp;
      std::sort(_alltc.begin(),_alltc.end(),comp);
    }
  }

  void TimeClusterDiag::findPrimaries (art::Event const& event,
      std::vector<spcount>& primaries) {
    _alltc.clear();
    primaries.clear();
    for(size_t ic = 0;ic< _tccol->size(); ++ ic){
      TimeCluster const& tc = _tccol->at(ic);
      spcount primary;
      findPrimary(event, tc, primary);
      primaries.push_back(primary);
    }
  }

  void TimeClusterDiag::findPrimary(art::Event const& event,
      TimeCluster const& tc, spcount& primary ){
    vector<spcount> sct;
    for (auto ich : tc._strawHitIdxs) {
      std::vector<StrawDigiIndex> shids;
      _chcol->fillStrawDigiIndices(ich,shids);
      for(auto shid : shids ) {
        StrawDigiMC const& mcdigi = _mcdigis->at(shid);
        StrawEnd itdc;
        art::Ptr<SimParticle> sp = mcdigi.earlyStrawGasStep()->simParticle();
        bool found(false);
        for(size_t isp=0;isp<sct.size();++isp){
          // count direct daughter/parent as part the same particle
          if(sct[isp]._spp == sp ){
            found = true;
            sct[isp].append(sp);
            break;
          }
        }
        if(!found)sct.push_back(sp);
      }
    }
    // sort by # of contributions
    sort(sct.begin(),sct.end(),spcountcomp());
    if(sct.size() > 0){
      primary = sct[0];
    }
  }

  void TimeClusterDiag::fillClusterInfo (TimeCluster const& tc,
      spcount const& primary, TimeClusterInfo& tcinfo) {
    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    const Calorimeter* calo = ch.get();
    // simple entries

    //tcinfo._nhits = tp._strawHitIdxs.size();

    tcinfo._nhits = tc.nStrawHits();
    tcinfo._time  = tc._t0._t0;
    tcinfo._terr  = tc._t0._t0err;
    tcinfo._pos   = tc._pos;

    // calo info if available
    if(tc._caloCluster.isNonnull()){
      tcinfo._ecalo = tc._caloCluster->energyDep();
      tcinfo._tcalo = _ttcalc.caloClusterTime(*tc._caloCluster,_pitch);
      tcinfo._dtcalo = _ttcalc.caloClusterTime(*tc._caloCluster,_pitch) - tc._t0._t0;
      // calculate the cluster position.  Currently the Z is in disk coordinates and must be translated, FIXME!
      XYZVectorF cog = XYZVectorF(calo->geomUtil().mu2eToTracker(calo->geomUtil().diskFFToMu2e(tc._caloCluster->diskID(),tc._caloCluster->cog3Vector())));
      tcinfo._cog = cog;
    }
    // mc info
    if(_mcdiag){
      tcinfo._prifrac = float(primary._count)/tc.nStrawHits();
      art::Ptr<SimParticle> const& sp = primary._spp;
      if(sp.isNonnull()){
        tcinfo._pripdg = sp->pdgId();
        tcinfo._priproc = sp->creationCode();
        if(sp->genParticle().isNonnull())
          tcinfo._prigen = sp->genParticle()->generatorId().id();
      }
    }

    // hit summary
    tcinfo._minhtime = 1.0e9;
    tcinfo._maxhtime = 0.0;
    tcinfo._ncehits = 0;
    // count CE hits
    for (auto const& idx : tc._strawHitIdxs) {
      ComboHit const& ch = _chcol->at(idx);
      float ht = _ttcalc.comboHitTime(ch,_pitch);
      tcinfo._maxhtime = std::max( ht , tcinfo._maxhtime);
      tcinfo._minhtime = std::min( ht , tcinfo._minhtime);
      if(_mcdiag){
        // mc truth info
        bool conversion = _shfcol[idx].hasAllProperties(_cesel);
        if(conversion)tcinfo._ncehits += ch.nStrawHits();
      }
    }
    // look for overlaps
    for(size_t ic = 0;ic< _tccol->size(); ++ ic){
      if((int)ic != tcinfo._tcindex) // don't count myself
        tcinfo._maxover = max(tcinfo._maxover,(Float_t)TrkUtilities::overlap(tc,_tccol->at(ic)));
    }
  }

  void TimeClusterDiag::createDiagnostics() {
    art::ServiceHandle<art::TFileService> tfs;
    // time peak diagnostics
    _tcdiag=tfs->make<TTree>("tcdiag","time cluster diagnostics");
    // event info should use the TrkAna struct, FIXME!!
    _tcdiag->Branch("iev",&_iev,"iev/I");
    _tcdiag->Branch("ntc",&_ntc,"ntc/I");
    _tcdiag->Branch("mcmidt0",&_mcmidt0,"mcmidt0/F");
    _tcdiag->Branch("besttc",&_besttc,TimeClusterInfo::leafnames().c_str());
    if(_mcdiag){
      _tcdiag->Branch("ceclust",&_ceclust,MCClusterInfo::leafnames().c_str());
    }
    if(_diag > 1) _tcdiag->Branch("tchinfo",&_tchinfo);
    if(_diag > 2) _tcdiag->Branch("alltc",&_alltc);
  }

  }  // end namespace mu2e

  using mu2e::TimeClusterDiag;
  DEFINE_ART_MODULE(TimeClusterDiag)
