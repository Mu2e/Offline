//
// track reco chain diagnostics.  This module also makes counts
// that can be used to estimate trigger rates.
//
// Original author D. Brown (LBNL) October 2016

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
#include "TrkDiag/inc/TrkMCTools.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
// data
#include "RecoDataProducts/inc/TimeClusterCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/HelixSeedCollection.hh"
#include "RecoDataProducts/inc/KalSeedCollection.hh"
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
#include <map>
#include <functional>
#include <algorithm>
#include <iostream>
using namespace std; 
using CLHEP::Hep3Vector;

namespace mu2e {

  typedef KalSeedCollection::const_iterator KSI;
  typedef HelixSeedCollection::const_iterator HSI;
  typedef TimeClusterCollection::const_iterator TCI;
  typedef art::Ptr<SimParticle> SPP;

  class TrkRecoDiag : public art::EDAnalyzer {
    public:
      explicit TrkRecoDiag(fhicl::ParameterSet const& pset);
      virtual ~TrkRecoDiag();
      virtual void beginRun(art::Run const& run) override;
      // This is called for each event.
      virtual void analyze(art::Event const& e) override;
    private:
// config parameters
      int _diag;
      bool _mcdiag;
      int _mcgen, _mcproc, _mcpdg; // targets for MC match
      TrkFitFlag _goodkfinal, _goodkseed, _goodhelix, _goodtc; // define a good track
     // event object tags
      art::InputTag _shTag;
      art::InputTag _shpTag;
      art::InputTag _shfTag;
      art::InputTag _hsTag;
      art::InputTag _ksTag;
      art::InputTag _kfTag;
      art::InputTag _tcTag;
      art::InputTag _mcdigisTag;
      art::InputTag _vdmcstepsTag;
      // cache of event objects
      const StrawHitCollection* _shcol;
      const StrawHitPositionCollection* _shpcol;
      const StrawHitFlagCollection* _shfcol;
      const HelixSeedCollection* _hscol;
      const KalSeedCollection* _kscol;
      const KalSeedCollection* _kfcol;
      const TimeClusterCollection* _tccol;
      const StrawDigiMCCollection* _mcdigis;
      const StepPointMCCollection* _vdmcsteps;
      // time offsets
      SimParticleTimeOffset _toff;
      // cache of BField at 0,0,0
      double _bz0;
      // helper functions
      bool findData(const art::Event& e);
      bool conversion(size_t index);
      // TTree and branch variables
      TTree *_trdiag;
      Int_t _iev;
      TrkFitFlag _kff, _ksf, _hsf; // fit flags
      HelixVal _kfh, _ksh; // helces at first segment
      Float_t _kft0, _kft0err, _kst0, _kst0err, _hst0, _hst0err, _tct0, _tct0err; // t0
      Float_t _kfm, _ksm, _kfmerr, _ksmerr; // momentum
      Int_t _kfn, _kfna, _ksn, _ksna, _hsn, _hsna, _tcn; // hit counts
      RobustHelix _hsh;
      Float_t _mcmom;
      Int_t _nprimary, _pdg, _gen, _proc;

      // helper functions 
      void findMCParticle(SPP& spp);
      KSI findKalFinal();
      KSI findKalSeed(KSI const& ikf);
      HSI findHelixSeed(KSI const& iks);
      TCI findTimeCluster(HSI const& ihs);
      KSI findMCMatch(KalSeedCollection const& ksc);
      HSI findMCMatch(HelixSeedCollection const& hsc);
      TCI findMCMatch(TimeClusterCollection const& tcc);
      KSI findBestReco(KalSeedCollection const& ksc, TrkFitFlag const& goodreco);
      HSI findBestReco(HelixSeedCollection const& hsc, TrkFitFlag const& goodreco);
      TCI findBestReco(TimeClusterCollection const& tcc);
      void fillKalSeed(KSI const& iks);
      void fillKalFinal(KSI const& ikf);
      void fillHelixSeed(HSI const& ihs);
      void fillTimeCluster(TCI const& itc);
      void fillMCData(KSI const& ksi);
      void fillMCData(HSI const& hsi);
      void fillMCData(TCI const& tci);
  };

  TrkRecoDiag::~TrkRecoDiag() {
  }

  TrkRecoDiag::TrkRecoDiag(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _diag(pset.get<int>("DiagLevel",1)),
    _mcdiag(pset.get<bool>("MonteCarloDiag",true)),
    _mcgen(pset.get<int>("MCGenerator",2)),// default for conversion electron
    _mcproc(pset.get<int>("MCProcess",56)),
    _mcpdg(pset.get<int>("MCPDG",11)),
    _goodkfinal(pset.get<vector<string> >("GoodKalFinalFlag",vector<string>{"KalmanOK"})),
    _goodkseed(pset.get<vector<string> >("GoodKalSeedFlag",vector<string>{"SeedOK"})),
    _goodhelix(pset.get<vector<string> >("GoodHelixFlag",vector<string>{"HelixOK"})),
    _goodtc(pset.get<vector<string> >("GoodTimeClusterFlag",vector<string>{"HitsOK"})),
    _shTag(pset.get<string>("StrawHitCollectionTag","makeSH")),
    _shpTag(pset.get<string>("StrawHitPositionCollectionTag","MakeStereoHits")),
    _shfTag(pset.get<string>("StrawHitFlagCollectionTag","PosHelixFinder")),
    _hsTag(pset.get<string>("HelixSeedCollectionTag","PosHelixFinder")),
    _ksTag(pset.get<string>("KalSeedFitCollectionTag","KSFDeM")),
    _kfTag(pset.get<string>("KalFinalFitCollectionTag","KSFFDeM")),
    _mcdigisTag(pset.get<art::InputTag>("StrawDigiMCCollection","makeSH")),
    _vdmcstepsTag(pset.get<art::InputTag>("VDStepPointMCCollection","detectorFilter:virtualdetector")),
    _toff(pset.get<fhicl::ParameterSet>("TimeOffsets"))
  {
    if(_diag > 0){
      art::ServiceHandle<art::TFileService> tfs;
      _trdiag=tfs->make<TTree>("trdiag","Track Reconstruction Diagnostics");
      _trdiag->Branch("iev",&_iev,"iev/I");
      _trdiag->Branch("tct0",&_tct0,"tct0/F");
      _trdiag->Branch("tct0err",&_tct0err,"tct0err/F");
      _trdiag->Branch("tcn",&_tcn,"tcn/I");
      _trdiag->Branch("hsf",&_hsf);
      _trdiag->Branch("hsh",&_hsh);
      _trdiag->Branch("hst0",&_hst0,"hst0/F");
      _trdiag->Branch("hst0err",&_hst0err,"hst0err/F");
      _trdiag->Branch("hsn",&_hsn,"hsn/I");
      _trdiag->Branch("hsna",&_hsna,"hsna/I");
      _trdiag->Branch("ksf",&_ksf);
      _trdiag->Branch("ksh",&_ksh);
      _trdiag->Branch("kst0",&_kst0,"kst0/F");
      _trdiag->Branch("kst0err",&_kst0err,"kst0err/F");
      _trdiag->Branch("ksm",&_ksm,"ksm/F");
      _trdiag->Branch("ksmerr",&_ksmerr,"ksmerr/F");
      _trdiag->Branch("ksn",&_ksn,"ksn/I");
      _trdiag->Branch("ksna",&_ksna,"ksna/I");
      _trdiag->Branch("kff",&_kff);
      _trdiag->Branch("kfh",&_kfh);
      _trdiag->Branch("kft0",&_kft0,"kft0/F");
      _trdiag->Branch("kft0err",&_kft0err,"kft0err/F");
      _trdiag->Branch("kfm",&_kfm,"kfm/F");
      _trdiag->Branch("kfmerr",&_kfmerr,"kfmerr/F");
      _trdiag->Branch("kfn",&_kfn,"kfn/I");
      _trdiag->Branch("kfna",&_kfna,"kfna/I");
      if(_mcdiag){
      	_trdiag->Branch("mcmom",&_mcmom,"mcmom/F");
	_trdiag->Branch("nprimary",&_nprimary,"nprimary/I");
	_trdiag->Branch("pdg",&_pdg,"pdg/I");
	_trdiag->Branch("gen",&_gen,"gen/I");
	_trdiag->Branch("proc",&_proc,"proc/I");
      }
      if(_diag > 1){
	if(_mcdiag){
	}
      }
    }
  }

  void TrkRecoDiag::beginRun(art::Run const& run){
    // get bfield
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    Hep3Vector vpoint_mu2e = det->toMu2e(Hep3Vector(0.0,0.0,0.0));
    _bz0 = bfmgr->getBField(vpoint_mu2e).z();
  }

  void TrkRecoDiag::analyze(art::Event const& evt) {
    _iev=evt.id().event();
// find the data
    if(findData(evt)) {
// find the MC particle
      SPP bestpart;
      if(_mcdiag){
	findMCParticle(bestpart);
	_mcmom = bestpart->startMomentum().mag(); // should be at tracker entrance FIXME!
	_pdg = bestpart->pdgId();
	_proc = bestpart->realCreationCode();
	if(bestpart->genParticle().isNonnull() )
	  _gen = bestpart->genParticle()->generatorId().id();
      }
// start with the final fits.  Find the 'best' fit according to MC match or reco criteria
      auto ikfs = findKalFinal();
// fill info about this fit
      fillKalFinal(ikfs);
      if(_mcdiag && ikfs != _kfcol->end()) fillMCData(ikfs);
// find the matching seed fit (or the best seed fit if the final failed)
      auto ikss = findKalSeed(ikfs);
      fillKalSeed(ikss);
      if(_mcdiag && ikfs == _kfcol->end() && ikss != _kscol->end()) fillMCData(ikss);
  // follow down to the helix fit with the same logic
      auto ihs = findHelixSeed(ikss);
      fillHelixSeed(ihs);
      if(_mcdiag && ikfs == _kfcol->end() && ikss == _kscol->end() && ihs != _hscol->end()) fillMCData(ihs);
  // follow down to the time cluster with the same logic
      auto itc = findTimeCluster(ihs);
      fillTimeCluster(itc);
      if(_mcdiag && ikfs == _kfcol->end() && ikss == _kscol->end() && ihs == _hscol->end() && itc!= _tccol->end()) fillMCData(itc);
   // fill the tree
      _trdiag->Fill();
    } else
      cout << "TrkRecoDiag_module can't find data" << endl;
  }

  bool TrkRecoDiag::findData(const art::Event& evt){
    _shcol = 0;
    _shpcol = 0;
    _shfcol = 0;
    _tccol = 0;
    _hscol = 0;
    _kscol = 0;
    _kfcol = 0;
    _mcdigis = 0;
    _vdmcsteps = 0;
// nb: getValidHandle does the protection (exception) on handle validity so I don't have to
    auto shH = evt.getValidHandle<StrawHitCollection>(_shTag);
    _shcol = shH.product();
    auto shpH = evt.getValidHandle<StrawHitPositionCollection>(_shpTag);
    _shpcol = shpH.product();
    auto shfH = evt.getValidHandle<StrawHitFlagCollection>(_shfTag);
    _shfcol = shfH.product();
    auto tcH = evt.getValidHandle<TimeClusterCollection>(_tcTag);
    _tccol = tcH.product();
    auto hsH = evt.getValidHandle<HelixSeedCollection>(_hsTag);
    _hscol = hsH.product();
    auto ksH = evt.getValidHandle<KalSeedCollection>(_ksTag);
    _kscol = ksH.product();
    auto kfH = evt.getValidHandle<KalSeedCollection>(_kfTag);
    _kfcol = kfH.product();
    if(_mcdiag){
      auto mcdH = evt.getValidHandle<StrawDigiMCCollection>(_mcdigisTag);
      _mcdigis = mcdH.product();
      auto mcstepsH = evt.getValidHandle<StepPointMCCollection>(_vdmcstepsTag);
      _vdmcsteps = mcstepsH.product();
      // update time offsets
      _toff.updateMap(evt);
    }

    return _shcol != 0 && _shpcol != 0 && _shfcol != 0 && _hscol != 0 && _kscol!= 0 && _kfcol != 0
      && ((_mcdigis != 0 && _vdmcsteps != 0 ) || !_mcdiag);
  }


  KSI TrkRecoDiag::findKalFinal() {
    // assume failure
    KSI retval = _kfcol->end();
    // if we're using MC info, find the best match to the primary particle
    // if we're using MC info, find the best match to the primary particle
    if(_mcdiag)
      retval = findMCMatch(*_kscol);
    else 
      retval = findBestReco(*_kscol,_goodkfinal);
    return retval;
  }


  KSI TrkRecoDiag::findKalSeed(KSI const& ikf) {
    KSI retval = _kscol->end();
    // see if the input final fit is valid.  If so, find the maximum cluster overlap
    if(ikf != _kfcol->end()){
      double maxover(0.0);
      auto kf = *ikf; 
      for(auto iks = _kscol->begin(); iks != _kscol->end(); ++iks ) {
	double over = TrkUtilities::overlap(kf,*iks);
	if(over > maxover) {
	  maxover = over;
	  retval = iks;
	}
      }
    } else {
      // if we're using MC info, find the best match to the primary particle
      if(_mcdiag)
	retval = findMCMatch(*_kscol);
      else 
	retval = findBestReco(*_kscol,_goodkseed);
    } 
    return retval;
  }

  HSI TrkRecoDiag::findHelixSeed(KSI const& iks) {
    HSI retval = _hscol->end();
    // see if the input final fit is valid.  If so, find the maximum cluster overlap
    if(iks != _kscol->end()){
      double maxover(0.0);
      auto ks = *iks; 
      for(auto ihs = _hscol->begin(); ihs != _hscol->end(); ++ihs ) {
	double over = TrkUtilities::overlap(ks,*ihs);
	if(over > maxover) {
	  maxover = over;
	  retval = ihs;
	}
      }
    } else {
      // if we're using MC info, find the best match to the primary particle
      if(_mcdiag)
	retval = findMCMatch(*_hscol);
      else 
	retval = findBestReco(*_hscol,_goodhelix);
    } 
    return retval;
  }

  TCI TrkRecoDiag::findTimeCluster(HSI const& ihs) {
    TCI retval = _tccol->end();
    // see if the input final fit is valid.  If so, find the maximum cluster overlap
    if(ihs != _hscol->end()){
      double maxover(0.0);
      auto hs = *ihs; 
      for(auto itc = _tccol->begin(); itc != _tccol->end(); ++itc ) {
	double over = TrkUtilities::overlap(hs,*itc);
	if(over > maxover) {
	  maxover = over;
	  retval = itc;
	}
      }
    } else {
      // if we're using MC info, find the best match to the primary particle
      if(_mcdiag)
	retval = findMCMatch(*_tccol);
      else 
	retval = findBestReco(*_tccol);
    } 
    return retval;
  }

  KSI TrkRecoDiag::findMCMatch(KalSeedCollection const& ksc) {
    auto retval = ksc.end();
    unsigned nprimary = 0;
    SPP spp; 
    for(auto iks = ksc.begin(); iks != ksc.end(); ++iks) {
      vector<StrawHitIndex> hits;
      for(auto const& hhit : iks->hits())
	hits.push_back(hhit.index());
      unsigned nmc = TrkMCTools::primaryParticle(spp,hits,_mcdigis);
      if(spp.isNonnull() &&
	  spp->genParticle().isNonnull() &&
	  spp->genParticle()->generatorId().id() == _mcgen &&
	  spp->realCreationCode() == _mcproc &&
	  spp->pdgId() == _mcpdg &&
	  nmc > nprimary){
	retval = iks;
	nprimary = nmc;
      }
    }
    return retval;
  }

  HSI TrkRecoDiag::findMCMatch(HelixSeedCollection const& hsc) {
    auto retval = hsc.end();
    unsigned nprimary = 0;
    SPP spp; 
    for(auto ihs = hsc.begin(); ihs != hsc.end(); ++ihs) {
      vector<StrawHitIndex> hits;
      for(auto const& hhit : ihs->hits())
	hits.push_back(hhit.index());
      unsigned nmc = TrkMCTools::primaryParticle(spp,hits,_mcdigis);
      if(spp.isNonnull() &&
	  spp->genParticle().isNonnull() &&
	  spp->genParticle()->generatorId().id() == _mcgen &&
	  spp->realCreationCode() == _mcproc &&
	  spp->pdgId() == _mcpdg &&
	  nmc > nprimary){
	retval = ihs;
	nprimary = nmc;
      }
    }
    return retval;
  }

  TCI TrkRecoDiag::findMCMatch(TimeClusterCollection const& tcc) {
    auto retval = tcc.end();
    unsigned nprimary = 0;
    SPP spp; 
    for(auto itc = tcc.begin(); itc != tcc.end(); ++itc) {
      unsigned nmc = TrkMCTools::primaryParticle(spp,itc->hits(),_mcdigis);
      if(spp.isNonnull() &&
	  spp->genParticle().isNonnull() &&
	  spp->genParticle()->generatorId().id() == _mcgen &&
	  spp->realCreationCode() == _mcproc &&
	  spp->pdgId() == _mcpdg &&
	  nmc > nprimary){
	retval = itc;
	nprimary = nmc;
      }
    }
    return retval;
  }

  KSI TrkRecoDiag::findBestReco(KalSeedCollection const& ksc, TrkFitFlag const& goodreco) {
    auto retval = ksc.end();
    // take the fit with the highest momentum
    // should add qulaity info testing FIXME!
    double maxmom(0.0);
    for(auto iks = ksc.begin(); iks != ksc.end(); ++iks) {
      if(iks->status().hasAllProperties(goodreco) &&
	  iks->segments().front()._mom > maxmom){
	maxmom = iks->segments().front()._mom;
	retval = iks;
      }
    }
    return retval;
  }

  HSI TrkRecoDiag::findBestReco(HelixSeedCollection const& hsc, TrkFitFlag const& goodreco) {
    auto retval = hsc.end();
    // take the fit with the most active hits
    // should add qulaity info testing FIXME!
    unsigned maxnhits(0);
    for(auto ihs = hsc.begin(); ihs != hsc.end(); ++ihs) {
      unsigned nhits(0);
      for(auto hh : ihs->hits())
	if(!hh.flag().hasAnyProperty(StrawHitFlag::outlier)) ++nhits;
      if(ihs->status().hasAllProperties(goodreco) &&
	  nhits > maxnhits) {
	maxnhits = nhits;
	retval = ihs;
      }
    }
    return retval;
  }

  TCI TrkRecoDiag::findBestReco(TimeClusterCollection const& tcc) {
    auto retval = tcc.end();
    // take the fit with the most active hits
    // should add qulaity info testing FIXME!
    unsigned maxnhits(0);
    for(auto itc = tcc.begin(); itc != tcc.end(); ++itc) {
      if(itc->hits().size() >  maxnhits) {
	maxnhits = itc->hits().size();
	retval = itc;
      }
    }
    return retval;
  }

  void TrkRecoDiag::fillKalFinal(KSI const& ikf) {
    if(ikf != _kfcol->end()){
      auto kf = *ikf;
      // fill branches
      _kff = kf.status();
      _kft0 = kf.t0().t0();
      _kft0err = kf.t0().t0Err();
      KalSegment const& fseg = kf.segments().front();
      _kfh = fseg.helix();
      _kfm = fseg.mom();
      _kfmerr = fseg.momerr();
      _kfn = kf.hits().size();
      // count the active hits
      _kfna = 0;
      for(auto tsh : kf.hits())
	if(tsh.flag().hasAllProperties(StrawHitFlag::active))++_kfna;
    } else {
      // reset
      _kff = TrkFitFlag();
      _kfh = HelixVal();
      _kft0 = _kft0err = 0.0;
      _kfm = _kfmerr = 0.0;
      _kfn = _kfna = 0;
    }
  }

  void TrkRecoDiag::fillKalSeed(KSI const& iks) {
    if(iks != _kscol->end()){
      auto ks = *iks;
      // fill branches
      _ksf = ks.status();
      _kst0 = ks.t0().t0();
      _kst0err = ks.t0().t0Err();
      KalSegment const& fseg = ks.segments().front();
      _ksh = fseg.helix();
      _ksm = fseg.mom();
      _ksmerr = fseg.momerr();
      _ksn = ks.hits().size();
      // count the active hits
      _ksna = 0;
      for(auto tsh : ks.hits())
	if(tsh.flag().hasAllProperties(StrawHitFlag::active))++_ksna;
    } else {
      // reset
      _ksf = TrkFitFlag();
      _ksh = HelixVal();
      _kst0 = _kst0err = 0.0;
      _ksm = _ksmerr = 0.0;
      _ksn = _ksna = 0;
    }
  }

  void TrkRecoDiag::fillHelixSeed(HSI const& ihs) {
    if(ihs != _hscol->end()){
      auto hs = *ihs;
      // fill branches
      _hsf = hs.status();
      _hsh = hs.helix();
      _hst0 = hs.t0().t0();
      _hst0err = hs.t0().t0Err();
      _hsn = hs.hits().size();
      // count the active hits
      _hsna = 0;
      for(auto tsh : hs.hits())
	if(!tsh.flag().hasAnyProperty(StrawHitFlag::outlier))++_hsna;
    } else {
      // reset
      _hsf = TrkFitFlag();
      _hsh = RobustHelix();
      _hst0 = _hst0err = 0.0;
      _hsn = _hsna = 0;
    }
  }

  void TrkRecoDiag::fillTimeCluster(TCI const& itc) {
    if(itc != _tccol->end()){
      auto tc = *itc;
      // fill branches
      _tct0 = tc.t0().t0();
      _tct0err = tc.t0().t0Err();
      _tcn = tc.hits().size();
      // count the active hits
    } else {
      // reset
      _tct0 = _tct0err = 0.0;
      _tcn = 0;
    }
  }

  void TrkRecoDiag::findMCParticle(SPP& spp){
    // reset
    spp = SPP(); 
    std::map<SPP,unsigned> spmap; // count how many digis each good SimParticle has
    // loop over the digis and find the ones that match
    for(auto mcd : *_mcdigis) {
      SPP sp;
      TrkMCTools::simParticle(sp,mcd);
      if(sp.isNonnull() &&
	  sp->genParticle().isNonnull() &&
	  sp->genParticle()->generatorId().id() == _mcgen &&
	  sp->realCreationCode() == _mcproc &&
	  sp->pdgId() == _mcpdg){
	auto ifnd = spmap.find(sp);
	if(ifnd == spmap.end())
	  spmap[sp] = 1;// initialize
	else
	  ++(ifnd->second);
      }
    }
    // find particle with the highest count
    unsigned maxcount(0);
    for(auto isp =spmap.begin();isp != spmap.end(); ++isp) {
      if(isp->second > maxcount){
	maxcount = isp->second;
	spp = isp->first;
      }
    }
  }
// needs implementing FIXME!  
  void TrkRecoDiag::fillMCData(KSI const& ksi) {
  }
  void TrkRecoDiag::fillMCData(HSI const& hsi) {
  }
  void TrkRecoDiag::fillMCData(TCI const& tci) {
  }

}
using mu2e::TrkRecoDiag;
DEFINE_ART_MODULE(TrkRecoDiag);
