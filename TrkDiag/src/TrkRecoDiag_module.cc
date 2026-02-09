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
#include "art_root_io/TFileService.h"
// mu2e
#include "Offline/GeneralUtilities/inc/Angles.hh"
#include "Offline/Mu2eUtilities/inc/MVATools.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
// diagnostics
#include "Offline/TrkDiag/inc/TrkMCTools.hh"
#include "Offline/TrkReco/inc/TrkUtilities.hh"
// data
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/TrkQual.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/EventWeight.hh"
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
  typedef art::Ptr<KalSeed> KSP;
  typedef art::Ptr<HelixSeed> HSP;
  typedef art::Ptr<TimeCluster> TCP;
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
      TrkFitFlag _goodkf, _goodks, _goodhs; // define a good track
      // event object tags
      art::InputTag _hsTag;
      art::InputTag _ksTag;
      art::InputTag _kfTag;
      art::InputTag _tqTag;
      art::InputTag _tcTag;
      art::InputTag _chTag;
      art::InputTag _mcdigisTag;
      art::InputTag _primaryTag;
      art::InputTag _vdmcstepsTag;
      art::InputTag _beamWtModule;
      // cache of event objects
      const HelixSeedCollection* _hscol;
      const KalSeedCollection* _kscol;
      const KalSeedCollection* _kfcol;
      const TimeClusterCollection* _tccol;
      const ComboHitCollection* _chcol;
      const StrawDigiMCCollection* _mcdigis;
      const PrimaryParticle* _primary;
      const StepPointMCCollection* _vdmcsteps;
      const TrkQualCollection* _tqcol;
      // Virtual Detector IDs
      vector<int> _midvids, _entvids, _xitvids;
      // cache of BField at 0,0,0
      double _bz0;
      // helper functions
      bool findData(const art::Event& e);
      bool conversion(size_t index);
      // TTree and branch variables
      TTree *_trdiag;
      Int_t _iev;
      Float_t _beamwt;
      TrkFitFlag _kff, _ksf, _hsf; // fit flags
      HelixVal _kfh, _ksh; // helces at tracker entrance and exit
      Float_t _kft0, _kft0err, _kst0, _kst0err, _hst0, _hst0err, _tct0, _tct0err; // t0
      Float_t _kfm, _ksm, _kfmerr, _ksmerr; // momentum
      Float_t _kschisq, _kfchisq, _kscon, _kfcon; //chisquared and consistency
      Int_t _kfn, _kfna, _ksn, _ksna, _hsn, _hsna, _tcn; // hit counts
      Int_t _kfns; // straw counts
      Float_t _kftq; // track quality MVA
      Int_t _kfnap, _ksnap, _hsnap, _kfnp, _ksnp, _hsnp, _tcnp; // primary particle hit counts
      RobustHelix _hsh;
      Float_t _mcgenmom;
      Float_t _mcentmom, _mcmidmom, _mcxitmom;
      Float_t _mcentpz, _mcmidpz, _mcxitpz;
      Float_t _mcentt0, _mcmidt0, _mcxitt0;
      Int_t _nprimary, _pdg, _gen, _proc;
      UInt_t _ndigitot;
      UInt_t _ndigi, _nkfprimary, _nksprimary, _nhsprimary, _ntcprimary;
      RobustHelix _mch;

      // helper functions
      unsigned findMCParticle(SPP& spp);
      KSI findMCMatch(SPP const& spp,KalSeedCollection const& ksc,unsigned& nprimary);
      HSI findMCMatch(art::Event const& evt, SPP const& spp,HelixSeedCollection const& hsc,unsigned& nprimary);
      TCI findMCMatch(art::Event const& evt, SPP const& spp,TimeClusterCollection const& tcc,unsigned& nprimary);
      KSI findBestReco(KalSeedCollection const& ksc, TrkFitFlag const& goodreco);
      HSI findBestReco(HelixSeedCollection const& ksc, TrkFitFlag const& goodreco);
      TCI findBestReco(TimeClusterCollection const& tcc);
      void fillKalFinal(SPP const& spp,KSI const& kfi);
      void fillKalSeed(SPP const& spp,KalSeed const& ks);
      void fillHelixSeed(SPP const& spp,HelixSeed const& hs);
      void fillTimeCluster(SPP const& spp,TimeCluster const& tc);
      bool fillMCInfo(SPP const& pspp);
      void resetTTree();
  };

  TrkRecoDiag::~TrkRecoDiag() {
  }

  TrkRecoDiag::TrkRecoDiag(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _diag(pset.get<int>("DiagLevel",1)),
    _mcdiag(pset.get<bool>("MonteCarloDiag",true)),
    _goodkf(pset.get<vector<string> >("GoodKalFinalFlag",vector<string>{"KalmanOK"})),
    _goodks(pset.get<vector<string> >("GoodKalSeedFlag",vector<string>{"SeedOK"})),
    _goodhs(pset.get<vector<string> >("GoodHelixFlag",vector<string>{"HelixOK"})),
    _hsTag(pset.get<art::InputTag>("HelixSeedCollection","HelixFinder:Positive")),
    _ksTag(pset.get<art::InputTag>("KalSeedFitCollection")),
    _kfTag(pset.get<art::InputTag>("KalFinalFitCollection")),
    _tqTag(pset.get<art::InputTag>("TrkQualTag","TrkQualDeM")),
    _tcTag(pset.get<art::InputTag>("TimeClusterCollection","TimeClusterFinder")),
    _chTag(pset.get<art::InputTag>("ComboHitCollection","makePH")),
    _mcdigisTag(pset.get<art::InputTag>("StrawDigiMCCollection","makeSD")),
    _primaryTag(pset.get<art::InputTag>("PrimaryParticleTag","FindMCPrimary")),
    _vdmcstepsTag(pset.get<art::InputTag>("VDStepPointMCCollection","detectorFilter:virtualdetector")),
    _beamWtModule( pset.get<art::InputTag>("beamWeightModule","PBIWeight" ))
    {
      if(_diag > 0){
        art::ServiceHandle<art::TFileService> tfs;
        _trdiag=tfs->make<TTree>("trdiag","Track Reconstruction Diagnostics");
        _trdiag->Branch("iev",&_iev,"iev/I");
        _trdiag->Branch("tct0",&_tct0,"tct0/F");
        _trdiag->Branch("tct0err",&_tct0err,"tct0err/F");
        _trdiag->Branch("tcn",&_tcn,"tcn/I");
        _trdiag->Branch("hsf.",&_hsf);
        _trdiag->Branch("hsh.",&_hsh);
        _trdiag->Branch("hst0",&_hst0,"hst0/F");
        _trdiag->Branch("hst0err",&_hst0err,"hst0err/F");
        _trdiag->Branch("hsn",&_hsn,"hsn/I");
        _trdiag->Branch("hsna",&_hsna,"hsna/I");
        _trdiag->Branch("ksf.",&_ksf);
        _trdiag->Branch("ksh.",&_ksh);
        _trdiag->Branch("kst0",&_kst0,"kst0/F");
        _trdiag->Branch("kst0err",&_kst0err,"kst0err/F");
        _trdiag->Branch("ksm",&_ksm,"ksm/F");
        _trdiag->Branch("ksmerr",&_ksmerr,"ksmerr/F");
        _trdiag->Branch("kschisq",&_kschisq,"kschisq/F");
        _trdiag->Branch("kscon",&_kscon,"kscon/F");
        _trdiag->Branch("ksn",&_ksn,"ksn/I");
        _trdiag->Branch("ksna",&_ksna,"ksna/I");
        _trdiag->Branch("kff.",&_kff);
        _trdiag->Branch("kfh.",&_kfh);
        _trdiag->Branch("kft0",&_kft0,"kft0/F");
        _trdiag->Branch("kft0err",&_kft0err,"kft0err/F");
        _trdiag->Branch("kfm",&_kfm,"kfm/F");
        _trdiag->Branch("kfmerr",&_kfmerr,"kfmerr/F");
        _trdiag->Branch("kfchisq",&_kfchisq,"kfchisq/F");
        _trdiag->Branch("kfcon",&_kfcon,"kfcon/F");
        _trdiag->Branch("kftq",&_kftq,"kftq/F");
        _trdiag->Branch("kfn",&_kfn,"kfn/I");
        _trdiag->Branch("kfna",&_kfna,"kfna/I");
        _trdiag->Branch("kfns",&_kfns,"kfns/I");
        if(_mcdiag){
          _trdiag->Branch("beamwt",&_beamwt,"beamwt/F");
          _trdiag->Branch("ndigitot",&_ndigitot,"ndigitot/i");
          _trdiag->Branch("mcgenmom",&_mcgenmom,"mcgenmom/F");
          _trdiag->Branch("mcentmom",&_mcentmom,"mcentmom/F");
          _trdiag->Branch("mcentpz",&_mcentpz,"mcentpz/F");
          _trdiag->Branch("mcentt0",&_mcentt0,"mcentt0/F");
          _trdiag->Branch("mcmidmom",&_mcmidmom,"mcmidmom/F");
          _trdiag->Branch("mcmidpz",&_mcmidpz,"mcmidpz/F");
          _trdiag->Branch("mcmidt0",&_mcmidt0,"mcmidt0/F");
          _trdiag->Branch("mcxitmom",&_mcxitmom,"mcxitmom/F");
          _trdiag->Branch("mcxitpz",&_mcxitpz,"mcxitpz/F");
          _trdiag->Branch("mcxitt0",&_mcxitt0,"mcxitt0/F");
          _trdiag->Branch("mch.",&_mch);
          _trdiag->Branch("nprimary",&_nprimary,"nprimary/i");
          _trdiag->Branch("pdg",&_pdg,"pdg/I");
          _trdiag->Branch("gen",&_gen,"gen/I");
          _trdiag->Branch("proc",&_proc,"proc/I");
          _trdiag->Branch("ndigi",&_ndigi,"ndigi/i");
          _trdiag->Branch("nkfprimary",&_nkfprimary,"nkfprimary/i");
          _trdiag->Branch("nksprimary",&_nksprimary,"nksprimary/i");
          _trdiag->Branch("nhsprimary",&_nhsprimary,"nhsprimary/i");
          _trdiag->Branch("ntcprimary",&_ntcprimary,"ntcprimary/i");
          _trdiag->Branch("kfnp",&_kfnp,"kfnp/I");
          _trdiag->Branch("kfnap",&_kfnap,"kfnap/I");
          _trdiag->Branch("ksnp",&_ksnp,"ksnp/I");
          _trdiag->Branch("ksnap",&_ksnap,"ksnap/I");
          _trdiag->Branch("hsnp",&_hsnp,"hsnp/I");
          _trdiag->Branch("hsnap",&_hsnap,"hsnap/I");
          _trdiag->Branch("tcnp",&_tcnp,"tcnp/I");
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
    // multiple VIDs for the tracker midplane: can we ever fix this??
    _midvids.push_back(VirtualDetectorId::TT_Mid);
    _midvids.push_back(VirtualDetectorId::TT_MidInner);
    _entvids.push_back(VirtualDetectorId::TT_FrontHollow);
    _entvids.push_back(VirtualDetectorId::TT_FrontPA);
    _xitvids.push_back(VirtualDetectorId::TT_Back);
  }

  void TrkRecoDiag::analyze(art::Event const& evt) {
    _iev=evt.id().event();
    resetTTree();
    // find the data
    if(findData(evt)) {
      // find the MC particle
      SPP bestpart;
      if(_mcdiag){
        // basic event information
        _ndigitot = _mcdigis->size();
        art::Handle<EventWeight> beamWtHandle;
        evt.getByLabel(_beamWtModule, beamWtHandle);
        if(beamWtHandle.isValid())
          _beamwt = beamWtHandle->weight();
        // find the MC true particle that best meets the specified requirements
        _nprimary = findMCParticle(bestpart);
        if(bestpart.isNonnull()){
          _mcgenmom = bestpart->startMomentum().mag(); // momentum at production
          _pdg = bestpart->pdgId();
          _proc = bestpart->originParticle().creationCode();
          if(bestpart->genParticle().isNonnull() )
            _gen = bestpart->genParticle()->generatorId().id();
          _ndigi = TrkMCTools::countDigis(bestpart,_mcdigis);
          fillMCInfo(bestpart);
          // find the best reco matches to this particle.
          auto bestkf = findMCMatch(bestpart,*_kfcol, _nkfprimary);
          if(bestkf != _kfcol->end()){
            fillKalFinal(bestpart,bestkf);
          } else {
            auto bestks = findMCMatch(bestpart,*_kscol, _nksprimary);
            if(bestks != _kscol->end()){
              fillKalSeed(bestpart,*bestks);
            } else {
              auto besths = findMCMatch(evt, bestpart,*_hscol, _nhsprimary);
              if(besths != _hscol->end()){
                fillHelixSeed(bestpart,*besths);
              } else {
                auto besttc = findMCMatch(evt, bestpart,*_tccol, _ntcprimary);
                if(besttc != _tccol->end()){
                  fillTimeCluster(bestpart,*besttc);
                }
              }
            }
          }
        }
      }
      // if there's no MC match: see if a bkg fit is found
      if(bestpart.isNull()){
        auto ikf = findBestReco(*_kfcol,_goodkf);
        if(ikf != _kfcol->end()){
          fillKalFinal(bestpart,ikf);
        } else {
          auto iks = findBestReco(*_kscol,_goodks);
          if(iks != _kscol->end()){
            fillKalSeed(bestpart,*iks);
          } else {
            auto ihs = findBestReco(*_hscol,_goodhs);
            if(ihs != _hscol->end()){
              fillHelixSeed(bestpart,*ihs);
            } else {
              auto itc = findBestReco(*_tccol);
              if(itc != _tccol->end()){
                fillTimeCluster(bestpart,*itc);
              }
            }
          }
        }
      }
      // fill the tree
      _trdiag->Fill();
    } else
      std::cout << "TrkRecoDiag_module can't find data" << std::endl;
  }

  bool TrkRecoDiag::findData(const art::Event& evt){
    _tccol = 0;
    _hscol = 0;
    _chcol = 0;
    _kscol = 0;
    _kfcol = 0;
    _mcdigis = 0;
    _primary = 0;
    _vdmcsteps = 0;
    // nb: getValidHandle does the protection (exception) on handle validity so I don't have to
    auto tcH = evt.getValidHandle<TimeClusterCollection>(_tcTag);
    _tccol = tcH.product();
    auto chH = evt.getValidHandle<ComboHitCollection>(_chTag);
    _chcol = chH.product();
    auto tqH = evt.getValidHandle<TrkQualCollection>(_tqTag);
    _tqcol = tqH.product();
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
      auto pph = evt.getValidHandle<PrimaryParticle>(_primaryTag);
      _primary = pph.product();
    }
    return _hscol != 0 && _kscol!= 0 && _kfcol != 0 && _tccol != 0 && _tqcol != 0
      && ((_mcdigis != 0 && _vdmcsteps != 0 && _primary != 0 ) || !_mcdiag);
  }

  KSI TrkRecoDiag::findMCMatch(SPP const& bestspp,KalSeedCollection const& ksc,unsigned& nprimary) {
    auto retval = ksc.end();
    nprimary = 0;
    SPP spp;
    for(auto iks = ksc.begin(); iks != ksc.end(); ++iks) {
      vector<StrawHitIndex> hits;
      for(auto const& hhit : iks->hits())
        hits.push_back(hhit.index());
      unsigned nmc = TrkMCTools::primaryParticle(spp,hits,_mcdigis);
      if(spp == bestspp && nmc > nprimary){
        retval = iks;
        nprimary = nmc;
      }
    }
    return retval;
  }

  HSI TrkRecoDiag::findMCMatch(art::Event const&  evt, SPP const& bestspp,HelixSeedCollection const& hsc,unsigned& nprimary) {
    auto retval = hsc.end();
    nprimary = 0;
    SPP spp;
    for(auto ihs = hsc.begin(); ihs != hsc.end(); ++ihs) {
      std::vector<StrawDigiIndex> sdis;
      for(size_t ihit=0;ihit < ihs->hits().size();ihit++)
        ihs->hits().fillStrawDigiIndices(ihit,sdis);
      unsigned nmc = TrkMCTools::primaryParticle(spp,sdis,_mcdigis);
      if(spp == bestspp && nmc > nprimary){
        retval = ihs;
        nprimary = nmc;
      }
    }
    return retval;
  }

  TCI TrkRecoDiag::findMCMatch(art::Event const&  evt, SPP const& bestspp,TimeClusterCollection const& tcc,unsigned& nprimary) {
    auto retval = tcc.end();
    nprimary = 0;
    SPP spp;
    for(auto itc = tcc.begin(); itc != tcc.end(); ++itc) {
      // translate from ComboHit to StrawDigi indices
      std::vector<StrawDigiIndex> sdis;
      for(auto ihit : itc->hits())
        _chcol->fillStrawDigiIndices(ihit,sdis);
      unsigned nmc = TrkMCTools::primaryParticle(spp,sdis,_mcdigis);
      if(spp == bestspp && nmc > nprimary){
        retval = itc;
        nprimary = nmc;
      }
    }
    return retval;
  }

  KSI TrkRecoDiag::findBestReco(KalSeedCollection const& ksc, TrkFitFlag const& goodreco) {
    auto retval = ksc.end();
    // take the fit with the highest momentum
    // should add quality info testing FIXME!
    double maxmom(0.0);
    for(auto iks = ksc.begin(); iks != ksc.end(); ++iks) {
      double mom = iks->segments().front().mom();
      if(iks->status().hasAllProperties(goodreco) && mom > maxmom){
        maxmom = mom;
        retval = iks;
      }
    }
    return retval;
  }

  HSI TrkRecoDiag::findBestReco(HelixSeedCollection const& hsc, TrkFitFlag const& goodreco) {
    auto retval = hsc.end();
    // take the fit with the highest product momentum*nhits
    double maxqual(0.0);
    for(auto ihs = hsc.begin(); ihs != hsc.end(); ++ihs) {
      double qual = ihs->hits().nStrawHits()*abs(ihs->helix().lambda()*ihs->helix().radius());
      if(ihs->status().hasAllProperties(goodreco) && qual > maxqual){
        maxqual = qual;
        retval = ihs;
      }
    }
    return retval;
  }

  TCI TrkRecoDiag::findBestReco(TimeClusterCollection const& tcc) {
    auto retval = tcc.end();
    // take the cluster with the most hits.  Should add quality later?
    unsigned maxnhits(0);
    for(auto itc = tcc.begin(); itc != tcc.end(); ++itc) {
      if(itc->nStrawHits() > maxnhits){
        maxnhits = itc->nStrawHits();
        retval = itc;
      }
    }
    return retval;
  }

  void TrkRecoDiag::resetTTree() {
    // mc truth
    _beamwt = 1.0;
    _mcgenmom = _mcentmom = _mcentpz = _mcentt0 = -1.0;
    _mcmidmom = _mcmidpz = _mcmidt0 = _mcxitmom = _mcxitpz = _mcxitt0 = -1.0;
    _pdg = _gen = _proc = -1;
    _ndigitot = _ndigi = _nkfprimary = _nksprimary = _nhsprimary = _ntcprimary = 0;
    _mch = RobustHelix();
    // reset
    _kff = TrkFitFlag();
    _kfh = HelixVal();
    _kft0 = _kft0err = 0.0;
    _kfchisq = _kfcon = _kftq -1.0;
    _kfm = _kfmerr = 0.0;
    _kfn = _kfna = _kfnp = _kfnap = _kfns = 0;
    // reset
    _ksf = TrkFitFlag();
    _ksh = HelixVal();
    _kst0 = _kst0err = 0.0;
    _kschisq = _kscon = -1.0;
    _ksm = _ksmerr = 0.0;
    _ksn = _ksna = _ksnp = _ksnap = 0;
    // reset
    _hsf = TrkFitFlag();
    _hsh = RobustHelix();
    _hst0 = _hst0err = 0.0;
    _hsn = _hsna = _hsnp = _hsnap = 0;
    // reset
    _tct0 = _tct0err = 0.0;
    _tcn = _tcnp = 0;
  }

  void TrkRecoDiag::fillKalFinal(SPP const& spp,KSI const& kfi) {
    KalSeed const& kf = *kfi;
    // find the corresponding TrkQual
    int index = std::distance(_kfcol->begin(), kfi );
    if(index < 0 || index >= int(_tqcol->size()))
      throw cet::exception("DIAG")<<"mu2e::TrkRecoDiag: FinalFit data inconsistent"<< endl;
    // fill branches
    _kff = kf.status();
    _kft0 = kf.t0().t0();
    _kft0err = kf.t0().t0Err();
    _kfchisq = kf.chisquared();
    _kfcon = kf.fitConsistency();
    _kftq = _tqcol->at(index).MVAOutput();
    KalSegment const& fseg = kf.segments().front();
    _kfh = fseg.helix();
    _kfm = fseg.mom();
    _kfmerr = fseg.momerr();
    _kfn = kf.hits().size();
    // count the active hits
    for(auto tsh : kf.hits()){
      if(tsh.flag().hasAllProperties(StrawHitFlag::active))++_kfna;
      if(spp.isNonnull()){
        StrawDigiMC const& mcdigi = _mcdigis->at(tsh.index());
        auto const& spmcp = mcdigi.earlyStrawGasStep();
        if ( spmcp->simParticle() == spp){
          ++_kfnp;
          if(!tsh.flag().hasAnyProperty(StrawHitFlag::outlier))++_kfnap;
        }
      }
    }
    // count straws
    _kfns = kf.straws().size();
  }

  void TrkRecoDiag::fillKalSeed(SPP const& spp,KalSeed const& ks) {
    // fill branches
    _ksf = ks.status();
    _kst0 = ks.t0().t0();
    _kst0err = ks.t0().t0Err();
    _kschisq = ks.chisquared();
    _kscon = ks.fitConsistency();
    KalSegment const& fseg = ks.segments().front();
    _ksh = fseg.helix();
    _ksm = fseg.mom();
    _ksmerr = fseg.momerr();
    _ksn = ks.hits().size();
    // count the active hits
    for(auto tsh : ks.hits()){
      if(tsh.flag().hasAllProperties(StrawHitFlag::active))++_ksna;
      if(spp.isNonnull()){
        StrawDigiMC const& mcdigi = _mcdigis->at(tsh.index());
        auto spmcp = mcdigi.earlyStrawGasStep();
        if ( spmcp->simParticle() == spp){
          ++_ksnp;
          if(!tsh.flag().hasAnyProperty(StrawHitFlag::outlier))++_ksnap;
        }
      }
    }
    // folllow down the reco chain
    //    auto hsP = ks.helix();
    //    if(hsP.isNonnull()){
    //      fillHelixSeed(spp,*hsP);
    //    } else
    //      std::cout << "No helix seed found in kal seed!" << std::endl;
  }

  void TrkRecoDiag::fillHelixSeed(SPP const& spp,HelixSeed const& hs) {
    // fill branches
    _hsf = hs.status();
    _hsh = hs.helix();
    _hst0 = hs.t0().t0();
    _hst0err = hs.t0().t0Err();
    // count the active hits
    for(auto hhit : hs.hits()){
      _hsn += hhit.nStrawHits();
      if(!hhit.flag().hasAnyProperty(StrawHitFlag::outlier))_hsna += hhit.nStrawHits();
      if(spp.isNonnull()){
        StrawDigiMC const& mcdigi = _mcdigis->at(hhit.index());
        auto spmcp = mcdigi.earlyStrawGasStep();
        if ( spmcp->simParticle() == spp){
          _hsnp += hhit.nStrawHits();
          if(!hhit.flag().hasAnyProperty(StrawHitFlag::outlier))_hsnap += hhit.nStrawHits();
        }
      }
    }
    // go down the reco chain
    auto tcp = hs.timeCluster();
    if(tcp.isNonnull()){
      fillTimeCluster(spp,*tcp);
    } else
      std::cout << "No time cluster found in helix seed!" << std::endl;
  }

  void TrkRecoDiag::fillTimeCluster(SPP const& spp,TimeCluster const& tc) {
    // fill branches
    _tct0 = tc.t0().t0();
    _tct0err = tc.t0().t0Err();
    _tcn = tc.nStrawHits();
    // count matching hits
    for(auto tchit : tc.hits()){
      StrawDigiMC const& mcdigi = _mcdigis->at(tchit);
      auto spmcp = mcdigi.earlyStrawGasStep();
      if ( spmcp->simParticle() == spp){
        ++_tcnp;
      }
    }
  }

  unsigned TrkRecoDiag::findMCParticle(SPP& spp){
    // reset
    spp = SPP();
    // pick the most related particle with the most hits
    unsigned nprimary(0);
    for(auto psp : _primary->primarySimParticles()){
      // loop over the digis and find the ones that match
      unsigned count(0);
      MCRelationship mcrel;
      SPP sp;
      for(auto mcd : *_mcdigis) {
        MCRelationship prel(psp,mcd.earlyStrawGasStep()->simParticle());
        if(prel == mcrel)
          count++;
        else if(prel > mcrel){
          mcrel = prel;
          count = 1;
          sp = mcd.earlyStrawGasStep()->simParticle();
        }
      }
      if(count > nprimary){
        nprimary = count;
        spp = sp;
      }
    }
    return nprimary;
  }

  bool TrkRecoDiag::fillMCInfo(art::Ptr<SimParticle> const& pspp) {
    bool retval(false);
    if(pspp.isNonnull() && _vdmcsteps != 0){
      GeomHandle<DetectorSystem> det;
      GlobalConstantsHandle<ParticleDataList> pdt;
      // find the earliest step associated with this particle passing the tracker midplane
      cet::map_vector_key trkid = pspp->id();
      auto jmc = _vdmcsteps->end();
      for(auto imc = _vdmcsteps->begin();imc != _vdmcsteps->end(); ++imc ) {
        // find matching steps
        if(  imc->trackId() == trkid &&
            (find(_midvids.begin(),_midvids.end(),imc->volumeId()) != _midvids.end()) ) {
          //    cout << "Found matching step " << *imc << endl;
          if(jmc == _vdmcsteps->end() || imc->time() < jmc->time())
            jmc = imc;
        }
      }
      if(jmc != _vdmcsteps->end()){
        // get momentum and position from this point
        _mcmidmom = jmc->momentum().mag();
        _mcmidpz = jmc->momentum().z();
        _mcmidt0 = jmc->time();
        Hep3Vector pos = det->toDetector(jmc->position());
        double charge = pdt->particle(pspp->pdgId()).charge();
        TrkUtilities::RobustHelixFromMom(pos,jmc->momentum(),charge,_bz0,_mch);
        retval = true;
      }
      // look for entrance and exit as well
      jmc = _vdmcsteps->end();
      for(auto imc = _vdmcsteps->begin();imc != _vdmcsteps->end(); ++imc ) {
        // find matching steps
        if(  imc->trackId() == trkid &&
            (find(_entvids.begin(),_entvids.end(),imc->volumeId()) != _entvids.end()) ) {
          //    cout << "Found matching step " << *imc << endl;
          if(jmc == _vdmcsteps->end() || imc->time() < jmc->time())
            jmc = imc;
        }
      }
      if(jmc != _vdmcsteps->end()){
        // get momentum and position from this point
        _mcentmom = jmc->momentum().mag();
        _mcentpz = jmc->momentum().z();
        _mcentt0 = jmc->time();
      }
      jmc = _vdmcsteps->end();
      for(auto imc = _vdmcsteps->begin();imc != _vdmcsteps->end(); ++imc ) {
        // find matching steps
        if(  imc->trackId() == trkid &&
            (find(_xitvids.begin(),_xitvids.end(),imc->volumeId()) != _xitvids.end()) ) {
          //    cout << "Found matching step " << *imc << endl;
          if(jmc == _vdmcsteps->end() || imc->time() < jmc->time())
            jmc = imc;
        }
      }
      if(jmc != _vdmcsteps->end()){
        // get momentum and position from this point
        _mcxitmom = jmc->momentum().mag();
        _mcxitpz = jmc->momentum().z();
        _mcxitt0 = jmc->time();
      }
    }
    return retval;
  }
}
using mu2e::TrkRecoDiag;
DEFINE_ART_MODULE(TrkRecoDiag)
