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
#include "MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "MCDataProducts/inc/EventWeight.hh"
#include "DataProducts/inc/threevec.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "TrkReco/inc/TrkUtilities.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"

// ROOT incldues
#include "Rtypes.h"

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
#include "TrkDiag/inc/KalDiag.hh"
#include "TrkDiag/inc/TrkComp.hh"
#include "TrkDiag/inc/HitCount.hh"
#include "TrkDiag/inc/TrkCount.hh"
#include "TrkDiag/inc/EventInfo.hh"
#include "TrkDiag/inc/TrkInfo.hh"
#include "TrkDiag/inc/EventWeightInfo.hh"
#include "TrkDiag/inc/TrkStrawHitInfo.hh"
#include "TrkDiag/inc/TrkStrawHitInfoMC.hh"
#include "TrkDiag/inc/TrkCaloHitInfo.hh"
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
    //TrkCaloMatchingParameters FitDirection and Track Particle
    TrkFitDirection _sdir;
    // CRV info
    std::string _crvCoincidenceModuleLabel;
    std::string _crvCoincidenceMCModuleLabel;
    // analysis options
    bool _fillmc, _pempty, _crv, _filltrkqual;
    int _diag;
    // analysis parameters
    double _minReflectTime; // minimum time for a track to reflect in the gradient
    double _mingood; // minimum fraction of momentum that classes a digi as good
    // Kalman fit diagnostics
    KalDiag _kdiag;
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
    art::InputTag _comboHitTag;
    art::InputTag _strawHitFlagTag;
    TrkCaloHitInfo _detch;
    std::vector<TrkStrawMatInfo> _detsm;
    // MC truth branches
    TrkInfoMC _demc, _uemc, _dmmc;
    art::InputTag _strawDigiMCTag;
    art::Handle<StrawDigiMCCollection> _strawDigiMCHandle;
    art::InputTag _simParticleTag;
    art::Handle<SimParticleCollection> _simParticleHandle;
    art::InputTag _vdStepPointMCTag;
    art::Handle<StepPointMCCollection> _vdStepPointMCHandle;
    // detailed MC truth for the signal candidate
    SimParticleTimeOffset _toff;
    TrkInfoMCStep _demcgen;
    TrkInfoMCStep _demcent, _demcmid, _demcxit;
    std::vector<TrkStrawHitInfoMC> _detshmc;
    // test trkqual variable branches
    TrkQualInfo _trkQualInfo;
    TrkQualTestInfo _trkqualTest;
    // helper functions
    void findMCData(const art::Event& event);
    void fillMCSteps(KalDiag::TRACKERPOS tpos, TrkFitDirection const& fdir, SimParticle::key_type id, TrkInfoMCStep& tmcs);
    void fillEventInfo(const art::Event& event);
    bool findBestTrack(KalSeedCollection const& kcol, KalSeed& kseed, TrkQualCollection const& tqcol, TrkQual& tqual);
    void resetBranches();
    void fillTrkInfoMC(art::Ptr<SimParticle> spp, const KalSeed& kseed,TrkInfoMC& trkinfomc) const;

    void fillTrkInfoMCStep(art::Ptr<SimParticle> spp, TrkInfoMCStep& trkinfomcstep) const;
    void fillTrkInfoMCStep(MCStepItr const& imcs, TrkInfoMCStep& trkinfomcstep) const;
    void fillTrkInfoMCStep(Hep3Vector const& mom, Hep3Vector const& pos, double charge, TrkInfoMCStep& trkinfomcstep) const;
    void fillMCInfo(const KalSeed& kseed, bool track_found);

    bool findUpstreamTrack(KalSeedCollection const& kcol,const KalSeed& dekseed, KalSeed& uekseed);
    bool findMuonTrack(KalSeedCollection const& kcol,const KalSeed& dekseed, KalSeed& dmukseed);
    // CRV info
    std::vector<CrvHitInfoReco> _crvinfo;
    std::vector<CrvHitInfoMC> _crvinfomc;
    // TestTrkQual
    void fillTrkQualInfo(const TrkQual& tqual, TrkQualInfo& trkqualInfo);
  };

  TrackAnalysisReco::TrackAnalysisReco(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset),
    _detag( pset.get<art::InputTag>("DeTag", art::InputTag()) ),
    _uetag( pset.get<art::InputTag>("UeTag", art::InputTag()) ),
    _dmtag( pset.get<art::InputTag>("DmuTag", art::InputTag()) ),
    _detqtag( pset.get<art::InputTag>("DeTrkQualTag", art::InputTag()) ),
    _meanPBItag( pset.get<art::InputTag>("MeanBeamIntensity",art::InputTag()) ),
    _PBIwtTag( pset.get<art::InputTag>("PBIWeightTag",art::InputTag()) ),
    _sdir((TrkFitDirection::FitDirection)(pset.get<int>("TrkFitDirection", TrkFitDirection::downstream))),
    _crvCoincidenceModuleLabel(pset.get<string>("CrvCoincidenceModuleLabel")),
    _crvCoincidenceMCModuleLabel(pset.get<string>("CrvCoincidenceMCModuleLabel")),
    _fillmc(pset.get<bool>("FillMCInfo",true)),
    _pempty(pset.get<bool>("ProcessEmptyEvents",true)),
    _crv(pset.get<bool>("AnalyzeCRV",false)),
    _filltrkqual(pset.get<bool>("fillTrkQualInfo",false)),
    _diag(pset.get<int>("diagLevel",1)),
    _minReflectTime(pset.get<double>("MinimumReflectionTime",20)), // nsec
    _mingood(pset.get<double>("MinimumGoodMomentumFraction",0.9)),
    _kdiag(pset.get<fhicl::ParameterSet>("KalDiag",fhicl::ParameterSet())),
    _trkana(0),
    _meanPBI(0.0),
    _comboHitTag(pset.get<art::InputTag>("ComboHitCollection", "")),
    _strawHitFlagTag(pset.get<art::InputTag>("StrawHitFlagCollection", "")),
    _strawDigiMCTag(pset.get<art::InputTag>("StrawDigiMCCollection", "")),
    _simParticleTag(pset.get<art::InputTag>("SimParticleCollection", "")),
    _vdStepPointMCTag(pset.get<art::InputTag>("VDStepPointMCCollection", "")),
    _toff(pset.get<fhicl::ParameterSet>("TimeOffsets"))
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
// optionally add detailed branches
    if(_diag > 1){
      _trkana->Branch("detsh",&_detsh);
      _trkana->Branch("detch",&_detch,TrkCaloHitInfo::leafnames().c_str());
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
    _toff.updateMap(event); //update time maps

    // Get handle to downstream electron track collection.  This also creates the final set of hit flags
    art::Handle<KalSeedCollection> deH;
    event.getByLabel(_detag,deH);
    // std::cout << _detag << std::endl; //teste
    KalSeedCollection const& deC = *deH;
    art::Handle<StrawHitFlagCollection> shfH;
    event.getByLabel(_strawHitFlagTag,shfH);
    StrawHitFlagCollection const& shfC = *shfH;
    art::Handle<ComboHitCollection> comboHitHandle;
    event.getByLabel(_comboHitTag, comboHitHandle);
    ComboHitCollection const& comboHits = *comboHitHandle;
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

    // reset
    resetBranches();

    // find the best track
    KalSeed dekseed, uekseed, dmukseed;
    bool deK = false, ueK = false, dmK = false; // bool to store whether this track was found
    TrkQual tqual;
    deK = findBestTrack(deC, dekseed, tqcol, tqual);
    if(deK != 0 || _pempty) {
      // find everything we want
      if (deK != 0) {
	ueK = findUpstreamTrack(ueC,dekseed, uekseed);
	dmK = findMuonTrack(dmC,dekseed, dmukseed);
      }
      if(_fillmc) { // get MC product collections
	findMCData(event);
      }

      // now fill out everything
      fillEventInfo(event);
      TrkTools::fillHitCount(shfC, _hcnt);
      if(deK != 0){
	TrkTools::fillTrkInfo(dekseed,_deti);
	_deti._trkqual = tqual.MVAOutput();
	if(_diag > 1){
	  TrkTools::fillHitInfo(dekseed, comboHits, _detsh); //TODO
	  TrkTools::fillMatInfo(dekseed, _detsm); //TODO
	}
	if(ueK != 0){
	  TrkTools::fillTrkInfo(uekseed,_ueti);
	}
	if(dmK != 0){
	  TrkTools::fillTrkInfo(dmukseed,_dmti);
	}
	if (dekseed.hasCaloCluster()) {
	  TrkTools::fillCaloHitInfo(dekseed, _detch); // TODO
	  _tcnt._ndec = 1; // only 1 possible calo hit at the moment
	}
	if (_filltrkqual) {
	  fillTrkQualInfo(tqual, _trkQualInfo);
	}
      }
      // fill mC info associated with this track
      if(_fillmc) fillMCInfo(dekseed, deK); 

      // fill CRV info
      if(_crv) CRVAnalysis::FillCrvHitInfoCollections(_crvCoincidenceModuleLabel, _crvCoincidenceMCModuleLabel, event, _crvinfo, _crvinfomc);

      // fill this row in the TTree
      _trkana->Fill();
    }
  }

  void TrackAnalysisReco::findMCData(const art::Event& event) {
    event.getByLabel(_strawDigiMCTag, _strawDigiMCHandle);
    event.getByLabel(_simParticleTag, _simParticleHandle);
    event.getByLabel(_vdStepPointMCTag, _vdStepPointMCHandle);
  }

  bool TrackAnalysisReco::findBestTrack(KalSeedCollection const& kcol, KalSeed& kseed, TrkQualCollection const& tqcol, TrkQual& tqual) {
    _tcnt._nde = kcol.size();
    // find the higest momentum track
    // TODO: should be finding the most "signal-like" track
    double max_momentum = -9999;
    bool track_found = false;
    for(size_t i_trk = 0; i_trk < kcol.size(); ++i_trk) {
      const KalSeed& i_kseed = kcol.at(i_trk);
      if (i_kseed.segments().begin()->mom() > max_momentum) {
	kseed = i_kseed; // currently takes the first track
	tqual = tqcol.at(i_trk);
	track_found = true;
      }
    }
    return track_found;
  }

  bool TrackAnalysisReco::findUpstreamTrack(KalSeedCollection const& kcol,const KalSeed& dekseed, KalSeed& uekseed) {
    _tcnt._nue = kcol.size();
    // loop over upstream tracks and pick the best one (closest to momentum) that's earlier than the downstream track
    bool track_found = false;
    double demom = dekseed.segments().begin()->mom();
    double closest_momentum = 0;
    for(auto i_kseed : kcol) {
      if(i_kseed.t0().t0() < dekseed.t0().t0() - _minReflectTime){
	double this_ue_momentum = i_kseed.segments().begin()->mom();
// choose the upstream track whose parameters best match the downstream track.
// Currently compare momentum at the tracker center, this should be done at the tracker entrance
// and should compare more parameters FIXME!
	if( fabs(this_ue_momentum-demom) < fabs(closest_momentum-demom)) {
	  uekseed = i_kseed;
	  track_found = true;
	}
      }
    }
    return track_found;
  }

  bool TrackAnalysisReco::findMuonTrack(KalSeedCollection const& kcol,const KalSeed& dekseed, KalSeed& dmukseed) {
    _tcnt._ndm = kcol.size();
    bool track_found = false;
// loop over muon tracks and pick the one with the largest hit overlap
    unsigned maxnover(0);
    for(auto i_kseed : kcol) {
      unsigned nover = _tcomp.nOverlap(i_kseed,dekseed);
      if(nover > maxnover){
	maxnover = nover;
	dmukseed = i_kseed;
	track_found = true;
      }
    }
    _tcnt._ndmo = maxnover;
    return track_found;
  }

  void TrackAnalysisReco::fillMCInfo(const KalSeed& kseed, bool track_found) {
  // for now MC truth is only for the DE track.  Maybe we need MC truth on other tracks too?  FIXME!
    art::Ptr<SimParticle> deSP;
    if (track_found) {
      TrkMCTools::findMCTrk(kseed, deSP, *_strawDigiMCHandle);
    } else if(_simParticleHandle.isValid()) { 
      // assume the 1st primary particle is the CE.  FIXME!!!
      const auto& simParticles = *_simParticleHandle;
      for ( auto isp = simParticles.begin(); isp != simParticles.end(); ++isp ){
    	if (isp->second.isSecondary()){
	  if (isp->second.parent()->isPrimary()){
	    deSP = isp->second.parent();
	    break;
	  }
	}
      }
    }
    if(deSP.isNonnull()){
      fillTrkInfoMC(deSP,kseed,_demc);
      fillTrkInfoMCStep(deSP,_demcgen);
      
      // find virtual detector steps where the particle crosses the tracker at fixed points
      static TrkFitDirection downstream;
      fillMCSteps(KalDiag::trackerEnt, downstream, cet::map_vector_key(deSP.key()), _demcent);
      fillMCSteps(KalDiag::trackerMid, downstream, cet::map_vector_key(deSP.key()), _demcmid);
      fillMCSteps(KalDiag::trackerExit, downstream, cet::map_vector_key(deSP.key()), _demcxit);
      if(_diag > 1 && track_found != 0) {
	// MC truth hit information
	TrkMCTools::fillHitInfoMCs(kseed, deSP, *_strawDigiMCHandle, _toff, _detshmc);
      }
    } 
  }

  void TrackAnalysisReco::fillMCSteps(KalDiag::TRACKERPOS tpos, TrkFitDirection const& fdir,
  SimParticle::key_type id, TrkInfoMCStep& tmcs) {
    std::vector<MCStepItr> steps;
    TrkMCTools::findMCSteps(*_vdStepPointMCHandle,id,_kdiag.VDids(tpos),steps);
    // if there are more than 1 crossing take the first one pointing in the specified direction
    for(auto istep = steps.begin(); istep != steps.end();++istep){
      if(std::signbit((*istep)->momentum().z()) == std::signbit(fdir.dzdt())){
	fillTrkInfoMCStep(*istep,tmcs);
	break;
      }
    }
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
    // clear vectors
    _detsh.clear();
    _detsm.clear();
    _detshmc.clear();
    _crvinfo.clear();
    _crvinfomc.clear();
  }


  void TrackAnalysisReco::fillTrkInfoMC(art::Ptr<SimParticle> spp, const KalSeed& kseed,TrkInfoMC& trkinfomc) const {
    // basic information
    if(spp->genParticle().isNonnull()) {
      trkinfomc._gen = spp->genParticle()->generatorId().id();
    }
    trkinfomc._pdg = spp->pdgId();
    trkinfomc._proc = spp->originParticle().creationCode();
    art::Ptr<SimParticle> pp = spp->originParticle().parent();
    if(pp.isNonnull()){
      trkinfomc._ppdg = pp->pdgId();
      trkinfomc._pproc = pp->originParticle().creationCode();
      trkinfomc._pmom = pp->startMomentum().vect().mag();
      if(pp->genParticle().isNonnull()) {
	trkinfomc._pgen = pp->genParticle()->generatorId().id();
      }
    }
    // fill track-specific  MC info
    int nactive = -1, nhits = -1, ngood = -1, nambig = -1;
    TrkMCTools::countHits(kseed, spp, *_strawDigiMCHandle, _mingood, nactive, nhits, ngood, nambig);
    trkinfomc._nactive = nactive;
    trkinfomc._nhits = nhits;
    trkinfomc._ngood = ngood;
    trkinfomc._nambig = nambig;


    int ndigi = -1, ndigigood = -1;
    TrkMCTools::countDigis(spp, *_strawDigiMCHandle, _mingood, ndigi, ndigigood);
    trkinfomc._ndigi = ndigi;
    trkinfomc._ndigigood = ndigigood;
  }

  void TrackAnalysisReco::fillTrkInfoMCStep(art::Ptr<SimParticle> spp, TrkInfoMCStep& trkinfomcstep) const {
    GlobalConstantsHandle<ParticleDataTable> pdt;
    GeomHandle<DetectorSystem> det;
    trkinfomcstep._time = _toff.totalTimeOffset(spp) + spp->startGlobalTime();
    double charge = pdt->particle(spp->pdgId()).ref().charge();
    Hep3Vector mom = spp->startMomentum();
    // need to transform into the tracker coordinate system
    Hep3Vector pos = det->toDetector(spp->startPosition());
    fillTrkInfoMCStep(mom,pos,charge,trkinfomcstep);
  }

  void TrackAnalysisReco::fillTrkInfoMCStep(MCStepItr const& imcs, TrkInfoMCStep& trkinfomcstep) const {
    GlobalConstantsHandle<ParticleDataTable> pdt;
    GeomHandle<DetectorSystem> det;
    trkinfomcstep._time = _toff.timeWithOffsetsApplied(*imcs);
    double charge = pdt->particle(imcs->simParticle()->pdgId()).ref().charge();
    Hep3Vector mom = imcs->momentum();
    // need to transform into the tracker coordinate system
    Hep3Vector pos = det->toDetector(imcs->position());
    fillTrkInfoMCStep(mom,pos,charge,trkinfomcstep);
  }

  void TrackAnalysisReco::fillTrkInfoMCStep(Hep3Vector const& mom, Hep3Vector const& pos, double charge, TrkInfoMCStep& trkinfomcstep) const {
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;

    trkinfomcstep._mom = mom.mag();
    trkinfomcstep._pos = pos;
    double hflt(0.0);
    HepVector parvec(5,0);
    static Hep3Vector vpoint_mu2e = det->toMu2e(Hep3Vector(0.0,0.0,0.0));
    static double bz = bfmgr->getBField(vpoint_mu2e).z();
    HepPoint ppos(pos.x(),pos.y(),pos.z());
    TrkHelixUtils::helixFromMom( parvec, hflt,ppos, mom,charge,bz);
    trkinfomcstep._hpar = helixpar(parvec);
  }


  void TrackAnalysisReco::fillTrkQualInfo(const TrkQual& tqual, TrkQualInfo& trkqualInfo) {
    int n_trkqual_vars = TrkQual::n_vars;
    for (int i_trkqual_var = 0; i_trkqual_var < n_trkqual_vars; ++i_trkqual_var) {
      TrkQual::MVA_varindex i_index = TrkQual::MVA_varindex(i_trkqual_var);
      trkqualInfo._trkqualvars[i_trkqual_var] = (double) tqual[i_index];
    }
    trkqualInfo._trkqual = tqual.MVAOutput();
  }
}  // end namespace mu2e

// Part of the magic that makes this class a module.
// create an instance of the module.  It also registers
using mu2e::TrackAnalysisReco;
DEFINE_ART_MODULE(TrackAnalysisReco);
