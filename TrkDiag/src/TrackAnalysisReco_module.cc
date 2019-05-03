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
#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "MCDataProducts/inc/EventWeight.hh"
#include "MCDataProducts/inc/KalSeedMC.hh"
#include "MCDataProducts/inc/CaloClusterMC.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "TrkReco/inc/TrkUtilities.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/TriggerResults.h"

// ROOT incldues
#include "Rtypes.h"
#include "TBits.h"
#include "TTree.h"
#include "TProfile.h"

// BaBar includes
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/TrkBase/TrkHelixUtils.hh"
#include "Mu2eUtilities/inc/TriggerResultsNavigator.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
// mu2e tracking
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
// diagnostics
#include "TrkDiag/inc/TrkComp.hh"
#include "TrkDiag/inc/HitCount.hh"
#include "TrkDiag/inc/TrkCount.hh"
#include "TrkDiag/inc/EventInfo.hh"
#include "TrkDiag/inc/TrkInfo.hh"
#include "TrkDiag/inc/GenInfo.hh"
#include "TrkDiag/inc/EventWeightInfo.hh"
#include "TrkDiag/inc/TrkStrawHitInfo.hh"
#include "TrkDiag/inc/TrkStrawHitInfoMC.hh"
#include "TrkDiag/inc/TrkCaloHitInfo.hh"
#include "TrkDiag/inc/CaloClusterInfoMC.hh"
#include "TrkDiag/inc/TrkQualInfo.hh"
#include "TrkDiag/inc/TrkQualTestInfo.hh"
#include "TrkDiag/inc/HelixInfo.hh"
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
    // reco count module
    art::InputTag _rctag;
    // CaloCrystal Ptr map
    art::InputTag _cchmtag;
    // SimParticleCollection Tag
    art::InputTag _spctag;
    //list of the triggerIds to track 
    std::vector<unsigned> _trigids;
    // event-weighting modules
    art::InputTag _meanPBItag;
    art::InputTag _PBIwtTag;
    // CRV info
    std::string _crvCoincidenceModuleLabel;
    std::string _crvCoincidenceMCModuleLabel;
    // analysis options
    bool _fillmc, _pempty, _crv, _helices, _filltrkqual, _filltrig;
    int _diag, _debug;
    // momentum analyzer
    double _bz0;
    // analysis parameters
    double _minReflectTime, _maxReflectTime; // minimum and maximum time for a track to reflect in the gradient
    // track comparator
    TrkComp _tcomp;
    // main TTree
    TTree* _trkana;
    TProfile* _tht; // profile plot of track hit times: just an example
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
    TrkCaloHitInfo _detch, _uetch;
    CaloClusterInfoMC _detchmc, _uetchmc;
    std::vector<TrkStrawMatInfo> _detsm;
    // trigger information
    unsigned _trigbits;
    // MC truth branches
    TrkInfoMC _demc, _uemc, _dmmc;
    art::InputTag _primaryParticleTag;
    art::InputTag _kalSeedMCTag, _caloClusterMCTag;
    std::vector<int> _entvids, _midvids, _xitvids;

    // detailed MC truth for the signal candidate
    GenInfo _demcgen, _demcpri; // generator and 'primary' information
    TrkInfoMCStep _demcent, _demcmid, _demcxit;
    std::vector<TrkStrawHitInfoMC> _detshmc;
    // test trkqual variable branches
    TrkQualInfo _trkQualInfo;
    TrkQualTestInfo _trkqualTest;
    // helper functions
    void fillEventInfo(const art::Event& event);
    void fillTriggerBits(const art::Event& event,std::string const& process);
//    TrkQualCollection const& tqcol, TrkQual& tqual);
    void resetBranches();
    KSCIter findBestRecoTrack(KalSeedCollection const& kcol);
    KSCIter findUpstreamTrack(KalSeedCollection const& kcol,KalSeed const& dekseed);
    KSCIter findMuonTrack(KalSeedCollection const& kcol,KalSeed const& dekseed);
    // CRV info
    std::vector<CrvHitInfoReco> _crvinfo;
    HelixInfo _hinfo;
    std::vector<CrvHitInfoMC> _crvinfomc;
    // SimParticle timing offset
    SimParticleTimeOffset _toff;
};

  TrackAnalysisReco::TrackAnalysisReco(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset),
    _detag( pset.get<art::InputTag>("DeTag", art::InputTag()) ),
    _uetag( pset.get<art::InputTag>("UeTag", art::InputTag()) ),
    _dmtag( pset.get<art::InputTag>("DmuTag", art::InputTag()) ),
    _detqtag( pset.get<art::InputTag>("DeTrkQualTag", art::InputTag()) ),
    _rctag( pset.get<art::InputTag>("RecoCountTag", art::InputTag()) ),
    _cchmtag( pset.get<art::InputTag>("CaloCrystalHitMapTag", art::InputTag()) ),
    _spctag( pset.get<art::InputTag>("SimParticleCollectionTag", art::InputTag()) ),
    _meanPBItag( pset.get<art::InputTag>("MeanBeamIntensity",art::InputTag()) ),
    _PBIwtTag( pset.get<art::InputTag>("PBIWeightTag",art::InputTag()) ),
    _crvCoincidenceModuleLabel(pset.get<string>("CrvCoincidenceModuleLabel")),
    _crvCoincidenceMCModuleLabel(pset.get<string>("CrvCoincidenceMCModuleLabel")),
    _fillmc(pset.get<bool>("FillMCInfo",true)),
    _pempty(pset.get<bool>("ProcessEmptyEvents",false)),
    _crv(pset.get<bool>("AnalyzeCRV",false)),
    _helices(pset.get<bool>("FillHelixInfo",false)),
    _filltrkqual(pset.get<bool>("FillTrkQualInfo",true)),
    _filltrig(pset.get<bool>("FillTriggerInfo",false)),
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _minReflectTime(pset.get<double>("MinimumReflectionTime",20)), // nsec
    _maxReflectTime(pset.get<double>("MaximumReflectionTime",200)), // nsec
    _trkana(0), _tht(0),
    _meanPBI(0.0),
    _primaryParticleTag(pset.get<art::InputTag>("PrimaryParticleTag", "")),
    _kalSeedMCTag(pset.get<art::InputTag>("KalSeedMCAssns", "")),
    _caloClusterMCTag(pset.get<art::InputTag>("CaloClusterMCAssns", "")),
    _toff(pset.get<fhicl::ParameterSet>("TimeOffsets"))
  {
    _midvids.push_back(VirtualDetectorId::TT_Mid);
    _midvids.push_back(VirtualDetectorId::TT_MidInner);
    _entvids.push_back(VirtualDetectorId::TT_FrontHollow);
    _entvids.push_back(VirtualDetectorId::TT_FrontPA);
    _xitvids.push_back(VirtualDetectorId::TT_Back);
  }

  void TrackAnalysisReco::beginJob( ){
    art::ServiceHandle<art::TFileService> tfs;
// create TTree
    _trkana=tfs->make<TTree>("trkana","track analysis");
    _tht=tfs->make<TProfile>("tht","Track Hit Time Profile",RecoCount::_nshtbins,-25.0,1725.0);
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
    if(_diag > 0){
      _trkana->Branch("uetch",&_uetch,TrkCaloHitInfo::leafnames().c_str());
    }
    if(_diag > 1){
      _trkana->Branch("detsh",&_detsh);
      _trkana->Branch("detsm",&_detsm);
    }
// add branches for other tracks
    _trkana->Branch("ue.",&_ueti,TrkInfo::leafnames().c_str());
    _trkana->Branch("dm.",&_dmti,TrkInfo::leafnames().c_str());
// trigger info.  Actual names should come from the BeginRun object FIXME
    if(_filltrig)_trkana->Branch("trigbits",&_trigbits,"trigbits/i");
// calorimeter information for the downstream electron track
// CRV info
   if(_crv) _trkana->Branch("crvinfo",&_crvinfo);
   // helix info
   if(_helices) _trkana->Branch("helixinfo",&_hinfo,HelixInfo::leafnames().c_str());
// optionally add MC truth branches
    if(_fillmc){
      _trkana->Branch("demc",&_demc,TrkInfoMC::leafnames().c_str());
      _trkana->Branch("demcgen",&_demcgen,GenInfo::leafnames().c_str());
      _trkana->Branch("demcpri",&_demcpri,GenInfo::leafnames().c_str());
      _trkana->Branch("demcent",&_demcent,TrkInfoMCStep::leafnames().c_str());
      _trkana->Branch("demcmid",&_demcmid,TrkInfoMCStep::leafnames().c_str());
      _trkana->Branch("demcxit",&_demcxit,TrkInfoMCStep::leafnames().c_str());
      if(_crv)_trkana->Branch("crvinfomc",&_crvinfomc);
      _trkana->Branch("detchmc",&_detchmc,CaloClusterInfoMC::leafnames().c_str());
      if(_diag > 0){
	_trkana->Branch("uetchmc",&_uetchmc,CaloClusterInfoMC::leafnames().c_str());
      }
      if(_diag > 1){
	_trkana->Branch("detshmc",&_detshmc);
      }
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
    // get bfield
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    Hep3Vector vpoint_mu2e = det->toMu2e(Hep3Vector(0.0,0.0,0.0));
    _bz0 = bfmgr->getBField(vpoint_mu2e).z();
  }

  void TrackAnalysisReco::analyze(const art::Event& event) {
  // update timing maps
    _toff.updateMap(event);
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
    // get the provenance from this for trigger processing
    std::string const& process = deH.provenance()->processName();
    // std::cout << _detag << std::endl; //teste
    auto const& deC = *deH;
    // find downstream muons and upstream electrons
    art::Handle<KalSeedCollection> ueH;
    event.getByLabel(_uetag,ueH);
    auto const& ueC = *ueH;
    art::Handle<KalSeedCollection> dmH;
    event.getByLabel(_dmtag,dmH);
    auto const& dmC = *dmH;
    art::Handle<CaloCrystalHitRemapping> cchmH;
    event.getByLabel(_cchmtag,cchmH);
    auto const& cchmap = *cchmH;
    art::Handle<SimParticleCollection> spcH;
    if(_fillmc){
      event.getByLabel(_spctag,spcH);
    }
    // general reco counts
    auto rch = event.getValidHandle<RecoCount>(_rctag);
    auto const& rc = *rch;
    for(size_t ibin=0;ibin < rc._nshtbins; ++ibin){
      float time = rc._shthist.binMid(ibin);
      float count  = rc._shthist.binContents(ibin);
      _tht->Fill(time,count);
    }
    // TrkQualCollection
    art::Handle<TrkQualCollection> trkQualHandle;
    event.getByLabel(_detqtag, trkQualHandle);
    TrkQualCollection const& tqcol = *trkQualHandle;
    // trigger information
    if(_filltrig){
      fillTriggerBits(event,process);
    }
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
    // process the best track
    if (idekseed != deC.end()) {
      auto const&  dekseed = *idekseed;
      TrkTools::fillTrkInfo(dekseed,_deti);
      if(_diag > 1){
	TrkTools::fillHitInfo(dekseed, _detsh);
	TrkTools::fillMatInfo(dekseed, _detsm);
      }
      if(_helices)TrkTools::fillHelixInfo(dekseed, _bz0, _hinfo);
      // upstream and muon tracks
      auto iuekseed = findUpstreamTrack(ueC,dekseed);
      if(iuekseed != ueC.end()) {
	auto const& uekseed = *iuekseed;
	TrkTools::fillTrkInfo(uekseed,_ueti);
	if(_diag >0 && uekseed.hasCaloCluster())
	  TrkTools::fillCaloHitInfo(uekseed, *caloh,  _uetch);
      }
	
      auto idmukseed = findMuonTrack(dmC,dekseed);
      if(idmukseed != dmC.end()) TrkTools::fillTrkInfo(*idmukseed,_dmti);
      // calorimeter info
      if (dekseed.hasCaloCluster()) {
	TrkTools::fillCaloHitInfo(dekseed, *caloh,  _detch);
	_tcnt._ndec = 1; // only 1 possible calo hit at the moment
	// test
	if(_debug>0){
	  auto const& tch = dekseed.caloHit();
	  auto const& cc = tch.caloCluster();
	  std::cout << "CaloCluster has energy " << cc->energyDep()
	  << " +- " << cc->energyDepErr() << std::endl;
	  for( auto const& cchptr: cc->caloCrystalHitsPtrVector() ) { 
	    // map the crystal ptr to the reduced collection
	    auto ifnd = cchmap.find(cchptr);
	    if(ifnd != cchmap.end()){
	      auto const& scchptr = ifnd->second;
	      if(scchptr.isNonnull())
		std::cout << "CaloCrystalHit has " << scchptr->energyDep() << " energy Dep" << std::endl;
	      else
		std::cout <<"CalCrystalHitPtr is invalid! "<< std::endl;
	    } else {
	      std::cout << "CaloCrystaLhitPtr not in map!" << std::endl;
	    }
	  }
	}
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
	//	std::cout << "KalSeedMCMatch has " << ksmcah->size() << " entries" << std::endl;
	for(auto iksmca = ksmcah->begin(); iksmca!= ksmcah->end(); iksmca++){
	//	  std::cout << "KalSeed Ptr " << dekptr << " match Ptr " << iksmca->first << std::endl;
	  if(iksmca->first == dekptr) {
	    auto const& dekseedmc = *(iksmca->second);
	    // primary associated SimParticle
	    auto trkprimary = dekseedmc.simParticle().simParticle(spcH);
	    TrkMCTools::fillTrkInfoMC(dekseedmc, trkprimary, dekseed, _demc);
	    double ttoff = _toff.totalTimeOffset(trkprimary); // kludge fix FIXME!
	    _demc._otime += ttoff; 
	    TrkMCTools::fillTrkInfoMCStep(dekseedmc, _demcent, _entvids);
	    TrkMCTools::fillTrkInfoMCStep(dekseedmc, _demcmid, _midvids);
	    TrkMCTools::fillTrkInfoMCStep(dekseedmc, _demcxit, _xitvids);
	    TrkMCTools::fillGenInfo(trkprimary, _demcgen, _demcpri, primary);
	    // times must be fixed FIXME!
	    _demcpri._time += ttoff;
	    _demcgen._time += ttoff;

	    if (_diag>1) {
	      TrkMCTools::fillHitInfoMCs(dekseedmc, _detshmc);
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
	if (_diag > 0 && iuekseed != ueC.end() && iuekseed->hasCaloCluster()) {
	  // fill MC truth of the associated CaloCluster 
	  for(auto iccmca= ccmcah->begin(); iccmca != ccmcah->end(); iccmca++){
	    if(iccmca->first == iuekseed->caloCluster()){
	      auto const& ccmc = *(iccmca->second);
	      TrkMCTools::fillCaloClusterInfoMC(ccmc,_uetchmc);
	      break;
	    }
	  }
	}
      }
    }
    if(idekseed != deC.end() || _pempty) {
      // fill general event information
      fillEventInfo(event);
      TrkTools::fillHitCount(rc, _hcnt);
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
      double this_momentum = kseed.segments().begin()->mom();
      if (this_momentum > max_momentum) {
	retval = i_kseed;
	max_momentum = this_momentum;
      }
    }
    return retval;
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

  void TrackAnalysisReco::fillTriggerBits(const art::Event& event,std::string const& process) {
    //get the TriggerResult from the process that created the KalFinalFit downstream collection
    art::InputTag const tag{Form("TriggerResults::%s", process.c_str())};
    auto trigResultsH = event.getValidHandle<art::TriggerResults>(tag);
    const art::TriggerResults* trigResults = trigResultsH.product();
    _trigbits = 0;
    for(size_t id=0;id < trigResults->size(); ++id){
      if (trigResults->accept(id)) {
	_trigbits |= 1 << id;
      }
    }
    if(_debug > 0){
      cout << "Found TriggerResults for process " << process << " with " << trigResults->size() << " Lines"
      << " trigger bits word " << _trigbits << endl;
      TriggerResultsNavigator tnav(trigResults);
      tnav.print();
    }
    
  }

  void TrackAnalysisReco::resetBranches() {
    // reset structs
    _einfo.reset();
    _hcnt.reset();
    _tcnt.reset();
    _deti.reset();
    _ueti.reset();
    _dmti.reset();
    _hinfo.reset();
    _demc.reset();
    _uemc.reset();
    _dmmc.reset();
    _demcgen.reset();
    _demcpri.reset();
    _demcent.reset();
    _demcmid.reset();
    _demcxit.reset();
    _wtinfo.reset();
    _trkqualTest.reset();
    _trkQualInfo.reset();
    _detch.reset();
    _uetch.reset();
    _detchmc.reset();
    _uetchmc.reset();
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
