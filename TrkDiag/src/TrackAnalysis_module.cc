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
//

// Mu2e includes
#include "GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "MCDataProducts/inc/EventWeight.hh"
#include "DataProducts/inc/threevec.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"

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
#include "BTrk/TrkBase/TrkParticle.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
// mu2e tracking
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
// diagnostics
#include "TrkDiag/inc/KalDiag.hh"
#include "TrkDiag/inc/TrkComp.hh"
#include "TrkDiag/inc/HitCount.hh"
#include "TrkDiag/inc/TrkCount.hh"
#include "TrkDiag/inc/TrkCaloDiag.hh"
#include "TrkDiag/inc/EventInfo.hh"
#include "TrkDiag/inc/TrkStrawHitInfo.hh"
#include "TrkDiag/inc/TrkStrawHitInfoMC.hh"
// CRV info
#include "CRVAnalysis/inc/CRVAnalysis.hh"

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
// This is fragile and needs to be last until CLHEP is
// properly qualified and included in the BaBar classes.
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"

using namespace std;

namespace mu2e {
// Need this for the BaBar headers.
  using CLHEP::Hep3Vector;

  class TrackAnalysis : public art::EDAnalyzer {

  public:
    explicit TrackAnalysis(fhicl::ParameterSet const& pset);
    virtual ~TrackAnalysis() { }

    void beginJob();
    void beginSubRun(const art::SubRun & subrun ) override;
    void analyze(const art::Event& e);

  private:

    // track collections.  Downstream electrons are signal candidates,
    // upstream electrons are used to identify cosmic background events
    // downstream muons are used in PID.  PID information should be analyzed
    // in a dedicated module FIXME!
    std::string _trkanaroot;
    art::InputTag _detag;
    art::InputTag _uetag;
    art::InputTag _dmtag;
    // event-weighting modules
    art::InputTag _genWttag;
    art::InputTag _PBItag, _meanPBItag;
    //TrkCaloMatchingParameters FitDirection and Track Particle
    TrkFitDirection _sdir;
    TrkParticle _spart;
    // CRV info
    std::string _crvCoincidenceModuleLabel;
    // analysis options
    bool _fillmc, _pempty, _crv;
    int _diag;
    // analysis parameters
    double _minReflectTime; // minimum time for a track to reflect in the gradient
    // Kalman fit diagnostics
    KalDiag _kdiag;
    // track comparator
    TrkComp _tcomp;
    // calorimeter diagnostics
    TrkCaloDiag _cdiag;
    TrkCaloInfo _dec;
    // main TTree
    TTree* _trkana;
    // general event info branch
    double _meanPBI;
    EventInfo _einfo;
    // hit counting
    HitCount _hcnt;
    // track counting
    TrkCount _tcnt;
    // track branches
    TrkInfo _deti, _ueti, _dmti;
    // detailed info branches for the signal candidate
    std::vector<TrkStrawHitInfo> _detsh;
    std::vector<TrkStrawMatInfo> _detsm;
    // MC truth branches
    TrkInfoMC _demc, _uemc, _dmmc;
    // detailed MC truth for the signal candidate
    TrkInfoMCStep _demcgen;
    TrkInfoMCStep _demcent, _demcmid, _demcxit;
    std::vector<TrkStrawHitInfoMC> _detshmc;
    // helper functions
    void fillMCSteps(KalDiag::TRACKERPOS tpos, TrkFitDirection const& fdir, SimParticle::key_type id, TrkInfoMCStep& tmcs);
    void fillEventInfo(const art::Event& event);
    void countHits(StrawHitFlagCollection const& shfC);
    const KalRep* findBestTrack(KalRepPtrCollection const& kcol);
    void resetBranches();
    void fillMCInfo(const KalRep* deK);
    void findBestClusterMatch(TrackClusterMatchCollection const& tcmc,
	const KalRep* krep, 
	TrackClusterMatchCollection::const_iterator& itcm);
    const KalRep* findUpstreamTrack(KalRepPtrCollection const& kcol,const KalRep* deK);
    const KalRep* findMuonTrack(KalRepPtrCollection const& kcol,const KalRep* deK);
    // CRV info
    std::vector<CrvHitInfoReco> _crvinfo;
    std::vector<CrvHitInfoMC> _crvinfomc;
  };
// instantiate statics

  TrackAnalysis::TrackAnalysis(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset),
    _trkanaroot(pset.get<std::string>("TrackAnalysisTagRoot","KFF") ),
//    _detag(pset.get<art::InputTag>("DownstreameTrack",art::InputTag()) ), // now determined from TrkParticle
//    _uetag(pset.get<art::InputTag>("UpstreameTrack",art::InputTag()) ), // now determined from TrkParticle
    _dmtag(pset.get<art::InputTag>("DownstreammuTrack",art::InputTag()) ),
    _genWttag( pset.get<art::InputTag>("generatorWeight",art::InputTag()) ),
    _PBItag( pset.get<art::InputTag>("BeamIntensity",art::InputTag()) ),
    _meanPBItag( pset.get<art::InputTag>("MeanBeamIntensity",art::InputTag()) ),
    _sdir((TrkFitDirection::FitDirection)(pset.get<int>("TrkFitDirection", TrkFitDirection::downstream))),
    _spart((TrkParticle::type)(pset.get<int>("TrkParticle"))),
    _crvCoincidenceModuleLabel(pset.get<string>("CrvCoincidenceModuleLabel")),
    _fillmc(pset.get<bool>("FillMCInfo",true)),
    _pempty(pset.get<bool>("ProcessEmptyEvents",true)),
    _crv(pset.get<bool>("AnalyzeCRV",false)),
    _diag(pset.get<int>("diagLevel",1)),
    _minReflectTime(pset.get<double>("MinimumReflectionTime",20)), // nsec
    _kdiag(pset.get<fhicl::ParameterSet>("KalDiag",fhicl::ParameterSet())),
    _cdiag(_spart,_sdir,pset.get<fhicl::ParameterSet>("TrkCaloDiag",fhicl::ParameterSet())),
    _trkana(0),
    _meanPBI(0.0)
  {
  }

  void TrackAnalysis::beginJob( ){
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
      _trkana->Branch("detsm",&_detsm);
    }
// add branches for other tracks
    _trkana->Branch("ue.",&_ueti,TrkInfo::leafnames().c_str());
    _trkana->Branch("dm.",&_dmti,TrkInfo::leafnames().c_str());
// calorimeter information for the downstream electron track
    _trkana->Branch("dec.",&_dec,TrkCaloInfo::leafnames().c_str());
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

  }

  void TrackAnalysis::beginSubRun(const art::SubRun & subrun ) {
    // mean number of protons on target
    art::Handle<ProtonBunchIntensity> PBIHandle;
    subrun.getByLabel(_meanPBItag, PBIHandle);
    if(PBIHandle.isValid())
      _meanPBI = PBIHandle->intensity();
  }
  void TrackAnalysis::analyze(const art::Event& event) {
    // Decide TrackTags from fit particle (muon tags specified in fcl)
    std::string chargename = _spart.charge() > 0.0 ? "P" : "M";
    _detag = _trkanaroot + "D" + _spart.name().substr(0,1) + chargename;
    _uetag = _trkanaroot + "U" + _spart.name().substr(0,1) + chargename;
    // Get handle to downstream electron track collection.  This also creates the final set of hit flags
    art::Handle<KalRepPtrCollection> deH;
    event.getByLabel(_detag,deH);
    // std::cout << _detag << std::endl; //teste
    KalRepPtrCollection const& deC = *deH;
    art::Handle<StrawHitFlagCollection> shfH;
    event.getByLabel(_detag,shfH);
    StrawHitFlagCollection const& shfC = *shfH;
    // find downstream muons and upstream electrons
    art::Handle<KalRepPtrCollection> ueH;
    event.getByLabel(_uetag,ueH);
    KalRepPtrCollection const& ueC = *ueH;
    art::Handle<KalRepPtrCollection> dmH;
    event.getByLabel(_dmtag,dmH);
    KalRepPtrCollection const& dmC = *dmH;
    // find Track-cluster matching data
    _cdiag.findData(event);
    // find the best track
    const KalRep* deK = findBestTrack(deC);
    if(deK != 0 || _pempty) {
      // reset
      resetBranches();
      // setup KalDiag.
      if(_fillmc)_kdiag.findMCData(event);
      // fill basic event information
      fillEventInfo(event);
      countHits(shfC);
      // fill the standard diagnostics
      if(deK != 0){
	_kdiag.fillTrkInfo(deK,_deti);
	if(_diag > 1){
	  _kdiag.fillHitInfo(deK, _detsh);
	  _kdiag.fillMatInfo(deK, _detsm);
	}
	// fill calorimeter information. First find the best matching cluster
	if(_cdiag.caloMatchHandle().isValid()){
	  TrackClusterMatchCollection const& tcmc = *_cdiag.caloMatchHandle();
	  TrackClusterMatchCollection::const_iterator itcm = tcmc.end();
	  findBestClusterMatch(tcmc,deK,itcm);
	  if(itcm != tcmc.end())
	    _cdiag.fillCaloInfo(*itcm,_dec);
	}
	// look for a matching upstream electron track
	const KalRep* ueK = findUpstreamTrack(ueC,deK);
	if(ueK != 0){
	  _kdiag.fillTrkInfo(ueK,_ueti);
	}
	// look for a matching muon track
	const KalRep* dmK = findMuonTrack(dmC,deK);
	if(dmK != 0){
	  _kdiag.fillTrkInfo(dmK,_dmti);
	}
      }
      // fill mC info associated with this track
      if(_fillmc) fillMCInfo(deK);

      // fill CRV info
      if(_crv) CRVAnalysis::FillCrvHitInfoCollections(_crvCoincidenceModuleLabel, event, _crvinfo, _crvinfomc);

      // fill this row in the TTree
      _trkana->Fill();
    }
  }

  const KalRep* TrackAnalysis::findBestTrack(KalRepPtrCollection const& kcol) {
    _tcnt._nde = kcol.size();
// if there aren't any tracks return 0
    const KalRep* deK(0);
    for(auto kptr : kcol ){
      const KalRep* krep = kptr.get();
      if(deK == 0) {
	deK = krep;
      } else {
// for now pick the highest-mometum track.
// This should pick out the most signal-like, best quailty track FIXME!!
	if(krep->momentum(0.0).mag() > deK->momentum(0.0).mag()){
// compute the hit overlap fraction
	  _tcnt._ndeo  = _tcomp.nOverlap(krep,deK);
	  deK = krep;
	}
      }
    }
    return deK;
  }

  const KalRep* TrackAnalysis::findUpstreamTrack(KalRepPtrCollection const& kcol,const KalRep* deK) {
    _tcnt._nue = kcol.size();
    const KalRep* ueK(0);
// loop over upstream tracks and pick the best one (closest to momentum) that's earlier than the downstream track
    for(auto kptr : kcol) {
      const KalRep* krep = kptr.get();
      if(krep->t0().t0() < deK->t0().t0() - _minReflectTime){
	if(ueK == 0){
	  ueK = krep;
	} else {
// choose the upstream track whose parameters best match the downstream track.
// Currently compare momentum at the tracker center, this should be done at the tracker entrance
// and should compare more parameters FIXME!
	  double demom = deK->momentum(0.0).mag();
	  if( fabs(krep->momentum(0.0).mag()-demom) < 
	      fabs(ueK->momentum(0.0).mag()-demom))
	    ueK = krep;
	}
      }
    }
    return ueK;
  }

  const KalRep* TrackAnalysis::findMuonTrack(KalRepPtrCollection const& kcol,const KalRep* deK) {
    _tcnt._ndm = kcol.size();
    const KalRep* dmK(0);
// loop over muon tracks and pick the one with the largest hit overlap
    unsigned maxnover(0);
    for(auto kptr : kcol) {
      const KalRep* krep = kptr.get();
      unsigned nover = _tcomp.nOverlap(krep,deK);
      if(nover > maxnover){
	maxnover = nover;
	dmK = krep;
      }
    }
    _tcnt._ndmo = maxnover;
    return dmK;
  }

  void TrackAnalysis::fillMCInfo(const KalRep* deK) {
  // for now MC truth is only for the DE track.  Maybe we need MC truth on other tracks too?  FIXME!
    art::Ptr<SimParticle> deSP;
    if(deK != 0) {
      _kdiag.findMCTrk(deK,deSP);
    } else if(_kdiag.mcData()._simparts != 0) { 
      // assume the 1st primary particle is the CE.  FIXME!!!
      for ( auto isp = _kdiag.mcData()._simparts->begin(); isp != _kdiag.mcData()._simparts->end(); ++isp ){
	/*	if(isp->second.isPrimary()){
	  deSP = art::Ptr<SimParticle>(_kdiag.mcData()._simparthandle,isp->second.id().asInt());
	  break;
	}
	*/
	if (isp->second.isSecondary()){
	  if (isp->second.parent()->isPrimary()){
	    deSP = isp->second.parent();
	    break;
	  }
	}
      }
    }
    if(deSP.isNonnull()){
      _kdiag.fillTrkInfoMC(deSP,deK,_demc);
      _kdiag.fillTrkInfoMCStep(deSP,_demcgen);
      // find virtual detector steps where the particle crosses the tracker at fixed points
      static TrkFitDirection downstream;
      fillMCSteps(KalDiag::trackerEnt, downstream, cet::map_vector_key(deSP.key()), _demcent);
      fillMCSteps(KalDiag::trackerMid, downstream, cet::map_vector_key(deSP.key()), _demcmid);
      fillMCSteps(KalDiag::trackerExit, downstream, cet::map_vector_key(deSP.key()), _demcxit);
      if(_diag > 1 && deK != 0) {
// MC truth hit information
	_kdiag.fillHitInfoMC(deSP, deK, _detshmc);
      }
    } 
  }

  void TrackAnalysis::fillMCSteps(KalDiag::TRACKERPOS tpos, TrkFitDirection const& fdir,
  SimParticle::key_type id, TrkInfoMCStep& tmcs) {
    std::vector<MCStepItr> steps;
    _kdiag.findMCSteps(_kdiag.mcData()._mcvdsteps,id,_kdiag.VDids(tpos),steps);
    // if there are more than 1 crossing take the first one pointing in the specified direction
    for(auto istep = steps.begin(); istep != steps.end();++istep){
      if(std::signbit((*istep)->momentum().z()) == std::signbit(fdir.dzdt())){
	_kdiag.fillTrkInfoMCStep(*istep,tmcs);
	break;
      }
    }
  }

  void TrackAnalysis::fillEventInfo( const art::Event& event) {
    // fill basic event information
    _einfo._eventid = event.event();
    _einfo._runid = event.run();
    _einfo._subrunid = event.subRun();
    // get event weight product
    _einfo._genwt = _einfo._beamwt = _einfo._evtwt = 1.; 
    _einfo._nprotons=-1;
    // total weight is the product of all weights
    // generator weight
    art::Handle<EventWeight> genWtHandle;
    event.getByLabel(_genWttag, genWtHandle);
    if(genWtHandle.isValid()){
      _einfo._genwt = genWtHandle->weight();
      _einfo._evtwt *= genWtHandle->weight();
    } 
    // actual number of protons on target
    art::Handle<ProtonBunchIntensity> PBIHandle;
    event.getByLabel(_PBItag, PBIHandle);
    if(PBIHandle.isValid()) {
      _einfo._nprotons = PBIHandle->intensity();
      if(_meanPBI > 0.0){
	_einfo._beamwt = _einfo._nprotons/_meanPBI; 
	_einfo._evtwt *= _einfo._beamwt;
      }
    }

  }

  void TrackAnalysis::findBestClusterMatch(TrackClusterMatchCollection const& tcmc,
  const KalRep* krep, 
  TrackClusterMatchCollection::const_iterator& itcm) {
    _tcnt._ndec = 0;
    itcm = tcmc.end();
    // for now pick highest-energy cluster.  This should match to the track or sum, FIXME!!
    double emin(-1.0);
    for( auto jtcm= tcmc.begin(); jtcm != tcmc.end(); ++jtcm ) {
      if(jtcm->textrapol()->trk().get() == krep){
	++_tcnt._ndec; // count number of matched clusters
	if(jtcm->caloCluster()->energyDep() > emin){
	  itcm = jtcm;
	  emin = jtcm->caloCluster()->energyDep(); 
	}
      }
    }
  }

  void TrackAnalysis::countHits(StrawHitFlagCollection const& shfC) {
    _hcnt._nsh = shfC.size();
    for(auto shf : shfC) {
      if(shf.hasAllProperties(StrawHitFlag::energysel))++_hcnt._nesel;
      if(shf.hasAllProperties(StrawHitFlag::radsel))++_hcnt._nrsel;
      if(shf.hasAllProperties(StrawHitFlag::timesel))++_hcnt._ntsel;
      if(shf.hasAllProperties(StrawHitFlag::bkg))++_hcnt._nbkg;
      if(shf.hasAllProperties(StrawHitFlag::stereo))++_hcnt._nster;
      if(shf.hasAllProperties(StrawHitFlag::tdiv))++_hcnt._ntdiv;
      if(shf.hasAllProperties(StrawHitFlag::trksel))++_hcnt._ntpk;
      if(shf.hasAllProperties(StrawHitFlag::elecxtalk))++_hcnt._nxt;
    }
  }

  void TrackAnalysis::resetBranches() {
  // reset structs
    _einfo.reset();
    _hcnt.reset();
    _tcnt.reset();
    _deti.reset();
    _ueti.reset();
    _dmti.reset();
    _dec.reset();
    _demc.reset();
    _uemc.reset();
    _dmmc.reset();
    _demcgen.reset();
    _demcent.reset();
    _demcmid.reset();
    _demcxit.reset();
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
using mu2e::TrackAnalysis;
DEFINE_ART_MODULE(TrackAnalysis);
