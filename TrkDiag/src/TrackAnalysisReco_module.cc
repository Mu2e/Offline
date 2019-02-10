//
// Prototype analysis module using tracks.  This module associates information from the
// Mu2e detector systems into a single coherent analysis TTree (trkana).  This module
// depends on the data products produced by reconstruction and (if requested) their MC
// counterparts.  The primary analysis object is the set of Downstream electron track fits.
// Upstream electron fit and downstream muon are also required for PID and quality selection.
// Calorimeter clusters and Track-cluster matching are used for PID. CRV coincidences are also
// included for rejecting cosmic backgrounds.
// Most of the calcluations are done by upstream modules and helper classes.
// Original author: Dave Brown (LBNL) 7/7/2016
// Updated November 2018 to run on KalSeeds only (A. Edmonds)
//

// Mu2e includes
#include "GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "MCDataProducts/inc/EventWeight.hh"
#include "MCDataProducts/inc/KalSeedMC.hh"
#include "MCDataProducts/inc/CaloClusterMC.hh"
#include "DataProducts/inc/threevec.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "TrkReco/inc/TrkUtilities.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"

// ROOT incldues
#include "Rtypes.h"
#include "TTree.h"

// BaBar includes
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/TrkBase/TrkHelixUtils.hh"
// mu2e tracking
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
// diagnostics
#include "TrkDiag/inc/TrkComp.hh"
#include "TrkDiag/inc/HitCount.hh"
#include "TrkDiag/inc/TrkCount.hh"
#include "TrkDiag/inc/EventInfo.hh"
#include "TrkDiag/inc/TrkInfo.hh"
#include "TrkDiag/inc/EventWeightInfo.hh"
#include "TrkDiag/inc/TrkStrawHitInfo.hh"
#include "TrkDiag/inc/TrkStrawHitInfoMC.hh"
#include "TrkDiag/inc/TrkCaloHitInfo.hh"
#include "TrkDiag/inc/CaloClusterInfoMC.hh"
#include "TrkDiag/inc/TrkQualInfo.hh"
#include "TrkDiag/inc/TrkQualTestInfo.hh"
#include "TrkDiag/inc/TrkTools.hh"
#include "TrkDiag/inc/TrkMCTools.hh"
// CRV info
#include "CRVAnalysis/inc/CRVAnalysis.hh"

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

using namespace std;

namespace mu2e {
// Need this for the BaBar headers.
  using CLHEP::Hep3Vector;
  typedef KalSeedCollection::const_iterator KSCIter;

  class TrackAnalysisReco : public art::EDAnalyzer {

  public:
    explicit TrackAnalysisReco(fhicl::ParameterSet const& pset);
    virtual ~TrackAnalysisReco() { }

    void beginJob();
    void beginSubRun(const art::SubRun & subrun ) override;
    void analyze(const art::Event& e);

  private:

    // track collections.  Downstream electrons are signal candidates,
    // upstream electrons are used to identify cosmic background events
    // downstream muons are used in PID.  PID information should be analyzed
    // in a dedicated module FIXME!
    art::InputTag _detag;
    art::InputTag _uetag;
    art::InputTag _dmtag;
    art::InputTag _detqtag;
    // event-weighting modules
    art::InputTag _meanPBItag;
    art::InputTag _PBIwtTag;
    // CRV info
    std::string _crvCoincidenceModuleLabel;
    std::string _crvCoincidenceMCModuleLabel;
    // analysis options
    bool _fillmc, _pempty, _crv, _filltrkqual;
    int _diag;
    // analysis parameters
    double _minReflectTime, _maxReflectTime; // minimum and maximum time for a track to reflect in the gradient
    // track comparator
    TrkComp _tcomp;
    // main TTree
    TTree* _trkana;
    // general event info branch
    double _meanPBI;
    EventInfo _einfo;
    EventWeightInfo _wtinfo;
    // hit counting
    HitCount _hcnt;
    // track counting
    TrkCount _tcnt;
    // track branches
    TrkInfo _deti, _ueti, _dmti;
    // detailed info branches for the signal candidate
    std::vector<TrkStrawHitInfo> _detsh;
    art::InputTag _strawHitFlagTag;
    TrkCaloHitInfo _detch;
    CaloClusterInfoMC _detchmc;
    std::vector<TrkStrawMatInfo> _detsm;
    // MC truth branches
    TrkInfoMC _demc, _uemc, _dmmc;
    art::InputTag _primaryParticleTag;
    art::InputTag _kalSeedMCTag, _caloClusterMCTag;
    // detailed MC truth for the signal candidate
    TrkInfoMCStep _demcgen;
    TrkInfoMCStep _demcent, _demcmid, _demcxit;
    std::vector<TrkStrawHitInfoMC> _detshmc;
    // test trkqual variable branches
    TrkQualInfo _trkQualInfo;
    TrkQualTestInfo _trkqualTest;
    // helper functions
    void fillEventInfo(const art::Event& event);
//    TrkQualCollection const& tqcol, TrkQual& tqual);
    void resetBranches();
    KSCIter findBestRecoTrack(KalSeedCollection const& kcol);
    KSCIter findPrimaryTrack(KalSeedCollection const& kcol,KalSeedMCAssns const& ksassn);
    KSCIter findUpstreamTrack(KalSeedCollection const& kcol,KalSeed const& dekseed);
    KSCIter findMuonTrack(KalSeedCollection const& kcol,KalSeed const& dekseed);
    // CRV info
    std::vector<CrvHitInfoReco> _crvinfo;
    std::vector<CrvHitInfoMC> _crvinfomc;
    // TestTrkQual
  };

  TrackAnalysisReco::TrackAnalysisReco(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset),
    _detag( pset.get<art::InputTag>("DeTag", art::InputTag()) ),
    _uetag( pset.get<art::InputTag>("UeTag", art::InputTag()) ),
    _dmtag( pset.get<art::InputTag>("DmuTag", art::InputTag()) ),
    _detqtag( pset.get<art::InputTag>("DeTrkQualTag", art::InputTag()) ),
    _meanPBItag( pset.get<art::InputTag>("MeanBeamIntensity",art::InputTag()) ),
    _PBIwtTag( pset.get<art::InputTag>("PBIWeightTag",art::InputTag()) ),
    _crvCoincidenceModuleLabel(pset.get<string>("CrvCoincidenceModuleLabel")),
    _crvCoincidenceMCModuleLabel(pset.get<string>("CrvCoincidenceMCModuleLabel")),
    _fillmc(pset.get<bool>("FillMCInfo",true)),
    _pempty(pset.get<bool>("ProcessEmptyEvents",true)),
    _crv(pset.get<bool>("AnalyzeCRV",false)),
    _filltrkqual(pset.get<bool>("fillTrkQualInfo",false)),
    _diag(pset.get<int>("diagLevel",1)),
    _minReflectTime(pset.get<double>("MinimumReflectionTime",20)), // nsec
    _maxReflectTime(pset.get<double>("MaximumReflectionTime",200)), // nsec
    _trkana(0),
    _meanPBI(0.0),
    _strawHitFlagTag(pset.get<art::InputTag>("StrawHitFlagCollection", "")),
    _primaryParticleTag(pset.get<art::InputTag>("PrimaryParticleTag", "")),
    _kalSeedMCTag(pset.get<art::InputTag>("KalSeedMCAssns", "")),
    _caloClusterMCTag(pset.get<art::InputTag>("CaloClusterMCAssns", ""))
  {
  }

  void TrackAnalysisReco::beginJob( ){
    art::ServiceHandle<art::TFileService> tfs;
// create TTree
    _trkana=tfs->make<TTree>("trkana","track analysis");
// add event info branch
    _trkana->Branch("evtinfo.",&_einfo,EventInfo::leafnames().c_str());
// hit counting branch
    _trkana->Branch("hcnt.",&_hcnt,HitCount::leafnames().c_str());
// track counting branch
    _trkana->Branch("tcnt.",&_tcnt,TrkCount::leafnames().c_str());
// add primary track (downstream electron) branch
    _trkana->Branch("de.",&_deti,TrkInfo::leafnames().c_str());
    //
    _trkana->Branch("detch",&_detch,TrkCaloHitInfo::leafnames().c_str());
// optionally add detailed branches
    if(_diag > 1){
      _trkana->Branch("detsh",&_detsh);
      _trkana->Branch("detsm",&_detsm);
    }
// add branches for other tracks
    _trkana->Branch("ue.",&_ueti,TrkInfo::leafnames().c_str());
    _trkana->Branch("dm.",&_dmti,TrkInfo::leafnames().c_str());
// calorimeter information for the downstream electron track
// CRV info
   if(_crv) _trkana->Branch("crvinfo",&_crvinfo);
// optionally add MC truth branches
    if(_fillmc){
      _trkana->Branch("demc",&_demc,TrkInfoMC::leafnames().c_str());
      _trkana->Branch("demcgen",&_demcgen,TrkInfoMCStep::leafnames().c_str());
      _trkana->Branch("demcent",&_demcent,TrkInfoMCStep::leafnames().c_str());
      _trkana->Branch("demcmid",&_demcmid,TrkInfoMCStep::leafnames().c_str());
      _trkana->Branch("demcxit",&_demcxit,TrkInfoMCStep::leafnames().c_str());
      if(_crv)_trkana->Branch("crvinfomc",&_crvinfomc);
      _trkana->Branch("detchmc",&_detchmc,CaloClusterInfoMC::leafnames().c_str());
      if(_diag > 1)_trkana->Branch("detshmc",&_detshmc);
    }
    if (_filltrkqual) {
      _trkana->Branch("detrkqual", &_trkQualInfo, TrkQualInfo::leafnames().c_str());
    }
  }

  void TrackAnalysisReco::beginSubRun(const art::SubRun & subrun ) {
    // mean number of protons on target
    art::Handle<ProtonBunchIntensity> PBIHandle;
    subrun.getByLabel(_meanPBItag, PBIHandle);
    if(PBIHandle.isValid())
      _meanPBI = PBIHandle->intensity();
  }

  void TrackAnalysisReco::analyze(const art::Event& event) {
  // get conditions/geometry objects
    mu2e::GeomHandle<mu2e::Calorimeter> caloh;
  // need to create and define the event weight branch here because we only now know the EventWeight creating modules that have been run through the Event
    if (!_trkana->GetBranch("evtwt")) { 
      std::vector<art::Handle<EventWeight> > eventWeightHandles;
      event.getManyByType(eventWeightHandles);
      if (eventWeightHandles.size()>0) {
	std::vector<std::string> labels;
	for (const auto& i_weightHandle : eventWeightHandles) {
	  std::string moduleLabel = i_weightHandle.provenance()->moduleLabel();
	  std::string instanceName = i_weightHandle.provenance()->productInstanceName();

	  std::string branchname = moduleLabel;
	  if (instanceName != "") {
	    branchname += "_" + instanceName;
	  }
	  labels.push_back(branchname);
	}
	_trkana->Branch("evtwt",&_wtinfo,_wtinfo.leafnames(labels).c_str());
      }
    }

    // Get handle to downstream electron track collection.  This also creates the final set of hit flags
    art::Handle<KalSeedCollection> deH;
    event.getByLabel(_detag,deH);
    // std::cout << _detag << std::endl; //teste
    KalSeedCollection const& deC = *deH;
    art::Handle<StrawHitFlagCollection> shfH;
    event.getByLabel(_strawHitFlagTag,shfH);
    StrawHitFlagCollection const& shfC = *shfH;
    // find downstream muons and upstream electrons
    art::Handle<KalSeedCollection> ueH;
    event.getByLabel(_uetag,ueH);
    KalSeedCollection const& ueC = *ueH;
    art::Handle<KalSeedCollection> dmH;
    event.getByLabel(_dmtag,dmH);
    KalSeedCollection const& dmC = *dmH;
    // TrkQualCollection
    art::Handle<TrkQualCollection> trkQualHandle;
    event.getByLabel(_detqtag, trkQualHandle);
    TrkQualCollection const& tqcol = *trkQualHandle;
// MC data
    art::Handle<PrimaryParticle> pph;
    art::Handle<KalSeedMCAssns> ksmcah;
    art::Handle<CaloClusterMCAssns> ccmcah;
    if(_fillmc) { // get MC product collections
      event.getByLabel(_primaryParticleTag,pph);
      event.getByLabel(_kalSeedMCTag,ksmcah);
      event.getByLabel(_caloClusterMCTag,ccmcah);
    }
    // reset
    resetBranches();
    // find the best tracks
    auto idekseed = findBestRecoTrack(deC);
    // if(_fillmc)idekseed = findPrimaryTrack(deC,*ksmcah);
    // process the best track
    if (idekseed != deC.end()) {
      auto const&  dekseed = *idekseed;
      TrkTools::fillTrkInfo(dekseed,_deti);
      if(_diag > 1){
	TrkTools::fillHitInfo(dekseed, _detsh); //TODO
	TrkTools::fillMatInfo(dekseed, _detsm); //TODO
      }
      // upstream and muon tracks
      auto iuekseed = findUpstreamTrack(ueC,dekseed);
      if(iuekseed != ueC.end()) TrkTools::fillTrkInfo(*iuekseed,_ueti);
      auto idmukseed = findMuonTrack(dmC,dekseed);
      if(idmukseed != dmC.end()) TrkTools::fillTrkInfo(*idmukseed,_dmti);
      // calorimeter info
      if (dekseed.hasCaloCluster()) {
	TrkTools::fillCaloHitInfo(dekseed, *caloh,  _detch); // TODO
	_tcnt._ndec = 1; // only 1 possible calo hit at the moment
      }
      if (_filltrkqual) {
	auto const& tqual = tqcol.at(std::distance(deC.begin(),idekseed));
	_deti._trkqual = tqual.MVAOutput();
	TrkTools::fillTrkQualInfo(tqual, _trkQualInfo);
      }
      // fill mC info associated with this track
      if(_fillmc ) { 
	const PrimaryParticle& primary = *pph;
	// use Assns interface to find the associated KalSeedMC; this uses ptrs
	auto dekptr = art::Ptr<KalSeed>(deH,std::distance(deC.begin(),idekseed));
	std::cout << "KalSeedMCMatch has " << ksmcah->size() << " entries" << std::endl;
	for(auto iksmca = ksmcah->begin(); iksmca!= ksmcah->end(); iksmca++){
	  std::cout << "KalSeed Ptr " << dekptr << " match Ptr " << iksmca->first << std::endl;
	  if(iksmca->first == dekptr) {
	    auto const& dekseedmc = *(iksmca->second);

	    TrkMCTools::fillTrkInfoMC(dekseedmc, dekseed, _demc);
	    TrkMCTools::fillTrkInfoMCStep(dekseedmc, _demcent, VirtualDetectorId::TT_FrontHollow); // TODO
	    TrkMCTools::fillTrkInfoMCStep(dekseedmc, _demcmid, VirtualDetectorId::TT_Mid); // TODO
	    TrkMCTools::fillTrkInfoMCStep(dekseedmc, _demcxit, VirtualDetectorId::TT_Back); // TODO

	    TrkMCTools::fillTrkInfoMCStep(dekseedmc, _demcgen, primary); // TODO
	    if (_diag>1) {
	      TrkMCTools::fillHitInfoMCs(dekseedmc, _detshmc); // TODO
	    }
	    break;
	  }
	}
	if (dekseed.hasCaloCluster()) {
	  // fill MC truth of the associated CaloCluster 
	  for(auto iccmca= ccmcah->begin(); iccmca != ccmcah->end(); iccmca++){
	    if(iccmca->first == dekseed.caloCluster()){
	      auto const& ccmc = *(iccmca->second);
	      TrkMCTools::fillCaloClusterInfoMC(ccmc,_detchmc);

	      break;
	    }
	  }
	}
      }
    }
    if(idekseed != deC.end() || _pempty) {
      // fill general event information
      fillEventInfo(event);
      TrkTools::fillHitCount(shfC, _hcnt);
      // TODO we want MC information when we don't have a track
      // fill CRV info
      if(_crv) CRVAnalysis::FillCrvHitInfoCollections(_crvCoincidenceModuleLabel, _crvCoincidenceMCModuleLabel, event, _crvinfo, _crvinfomc);
      // fill this row in the TTree
      _trkana->Fill();
    }
  }

  KSCIter TrackAnalysisReco::findBestRecoTrack(KalSeedCollection const& kcol) {
    KSCIter retval = kcol.end();
    _tcnt._nde = kcol.size();
    // find the higest momentum track; should be making some quality cuts too FIXME!
    double max_momentum = -9999;
    for(auto i_kseed=kcol.begin(); i_kseed != kcol.end(); ++i_kseed) {
      auto const& kseed = *i_kseed; 
      if (kseed.segments().begin()->mom() > max_momentum) {
	retval = i_kseed;
      }
    }
    return retval;
  }

  // find the track most associated with the MC true primary FIXME!
  KSCIter TrackAnalysisReco::findPrimaryTrack( KalSeedCollection const& kcol,
      KalSeedMCAssns const& mcassns) {
    return kcol.end();
  }

  KSCIter TrackAnalysisReco::findUpstreamTrack(KalSeedCollection const& kcol,const KalSeed& dekseed) {
    KSCIter retval = kcol.end();
    _tcnt._nue = kcol.size();
    // loop over upstream tracks and pick the best one (closest to momentum) that's earlier than the downstream track
    double demom = dekseed.segments().begin()->mom();
    double closest_momentum = 0;
    for(auto i_kseed=kcol.begin(); i_kseed != kcol.end(); i_kseed++) {
      if(i_kseed->t0().t0() < dekseed.t0().t0() - _minReflectTime &&
	i_kseed->t0().t0() > dekseed.t0().t0() - _maxReflectTime) {
	double this_ue_momentum = i_kseed->segments().begin()->mom();
// choose the upstream track whose parameters best match the downstream track.
// Currently compare momentum at the tracker center, this should be done at the tracker entrance
// and should compare more parameters FIXME!
	if( fabs(this_ue_momentum-demom) < fabs(closest_momentum-demom)) {
	  retval = i_kseed;
	}
      }
    }
    return retval;
  }

  KSCIter TrackAnalysisReco::findMuonTrack(KalSeedCollection const& kcol,const KalSeed& dekseed) {
    KSCIter retval = kcol.end();
    _tcnt._ndm = kcol.size();
// loop over muon tracks and pick the one with the largest hit overlap
    unsigned maxnover(0);
    for(auto i_kseed = kcol.begin(); i_kseed != kcol.end(); i_kseed++) {
      unsigned nover = _tcomp.nOverlap(*i_kseed,dekseed);
      if(nover > maxnover){
	maxnover = nover;
	retval = i_kseed;
      }
    }
    _tcnt._ndmo = maxnover;
    return retval;
  }

  void TrackAnalysisReco::fillEventInfo( const art::Event& event) {
    // fill basic event information
    _einfo._eventid = event.event();
    _einfo._runid = event.run();
    _einfo._subrunid = event.subRun();

    // get event weight products
    std::vector<art::Handle<EventWeight> > eventWeightHandles;
    event.getManyByType(eventWeightHandles);
    std::vector<Float_t> weights;
    for (const auto& i_weightHandle : eventWeightHandles) {
      double weight = i_weightHandle->weight();
      if (i_weightHandle.provenance()->moduleLabel() == _PBIwtTag.label()) {
	if (_meanPBI > 0.0){
	  _einfo._nprotons = _meanPBI*weight;
	}
	else {
	  _einfo._nprotons = 1; // for non-background mixed jobs
	}
      }
      weights.push_back(weight);
    }
    _wtinfo.setWeights(weights);
  }

  void TrackAnalysisReco::resetBranches() {
  // reset structs
    _einfo.reset();
    _hcnt.reset();
    _tcnt.reset();
    _deti.reset();
    _ueti.reset();
    _dmti.reset();
    _demc.reset();
    _uemc.reset();
    _dmmc.reset();
    _demcgen.reset();
    _demcent.reset();
    _demcmid.reset();
    _demcxit.reset();
    _wtinfo.reset();
    _trkqualTest.reset();
    _trkQualInfo.reset();
    _detch.reset();
    _detchmc.reset();
// clear vectors
    _detsh.clear();
    _detsm.clear();
    _detshmc.clear();
    _crvinfo.clear();
    _crvinfomc.clear();
  }

}  // end namespace mu2e

// Part of the magic that makes this class a module.
// create an instance of the module.  It also registers
using mu2e::TrackAnalysisReco;
DEFINE_ART_MODULE(TrackAnalysisReco);
