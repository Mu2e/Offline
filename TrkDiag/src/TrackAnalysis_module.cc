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
#include "TrkDiag/inc/KalDiag.hh"
#include "TrkDiag/inc/TrkCaloDiag.hh"
#include "TrkDiag/inc/EventInfo.hh"
#include "TrkDiag/inc/TrkHitShare.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
// C++ includes.
#include <iostream>
#include <string>

//G4Beamline includes
#include "MCDataProducts/inc/G4BeamlineInfoCollection.hh"

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
    void analyze(const art::Event& e);

  private:

    // track collections.  Downstream electrons are signal candidates,
    // upstream electrons and downstream muons are used to identify background events
    art::InputTag _demtag;
    art::InputTag _uemtag;
    art::InputTag _dmmtag;
    // event-weighting modules
    art::InputTag _genWttag;
    art::InputTag _beamWttag;
    art::InputTag _PBItag;
    vector<art::InputTag> _evtWttags;
    // Kalman fit diagnostics
    KalDiag _kdiag;
    // calorimeter diagnostics
    TrkCaloDiag _cdiag;
    TrkCaloInfo _demc;
    // options
    bool _fillmc, _diag;
    // main TTree
    TTree* _trkana;
    // general event info branch
    EventInfo _einfo;
    // track branches
    TrkInfo _demti, _uemti, _dmmti;
    // detailed info branches for the signal candidate
    std::vector<TrkStrawHitInfo> _demtsh;
    std::vector<TrkStrawMatInfo> _demtsm;
    // MC truth branches
    TrkInfoMC _demmc, _uemmc, _dmmmc;
    // detailed MC truth for the signal candidate
    TrkInfoMCStep _demmcgen;
    TrkInfoMCStep _demmcent;
    std::vector<TrkStrawHitInfoMC> _demtshmc;
    // helper functions
    void fillEventInfo(const art::Event& event);
    const KalRep* findBestTrack(KalRepPtrCollection const& kcol);
    void resetBranches();
    void fillMCInfo(const KalRep* demK);
    void findBestClusterMatch(TrackClusterMatchCollection const& tcmc,
	const KalRep* krep, 
	TrackClusterMatchCollection::const_iterator& itcm);
    // define signal direction and particle
    static TrkFitDirection _sdir;
    static TrkParticle _spart;
  };
// instantiate statics
  TrkFitDirection TrackAnalysis::_sdir(TrkFitDirection::downstream);
  TrkParticle TrackAnalysis::_spart(TrkParticle::e_minus);

  TrackAnalysis::TrackAnalysis(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset),
    _demtag(pset.get<art::InputTag>("DownstreameMinusTrackTag",art::InputTag()) ),
    _uemtag(pset.get<art::InputTag>("UpstreameMinusTrackTag",art::InputTag()) ),
    _dmmtag(pset.get<art::InputTag>("DownstreammuMinusTrackTag",art::InputTag()) ),
    _genWttag( pset.get<art::InputTag>("generatorWeightTag",art::InputTag()) ),
    _beamWttag( pset.get<art::InputTag>("beamWeightTag",art::InputTag()) ),
    _PBItag( pset.get<art::InputTag>("ProtonBunchIntensityTag",art::InputTag("ProtonBunchIntensitySummarizer")) ),
    _evtWttags( pset.get<std::vector<art::InputTag>>("eventWeightTags",std::vector<art::InputTag>() ) ),
    _kdiag(pset.get<fhicl::ParameterSet>("KalDiag",fhicl::ParameterSet())),
    _cdiag(_spart,_sdir,pset.get<fhicl::ParameterSet>("TrkCaloDiag",fhicl::ParameterSet())),
    _fillmc(pset.get<bool>("FillMCInfo",true)),
    _diag(pset.get<int>("diagLevel",1)),
    _trkana(0)
  {
  }

  void TrackAnalysis::beginJob( ){
    art::ServiceHandle<art::TFileService> tfs;
// create TTree
    _trkana=tfs->make<TTree>("trkana","track analysis");
// add event info branch
    _trkana->Branch("evtinfo",&_einfo,EventInfo::leafnames().c_str());
// add primary track (downstream electron) branch
    _trkana->Branch("dem",&_demti,TrkInfo::leafnames().c_str());
// optionally add detailed branches
    if(_diag > 1){
      _trkana->Branch("demtsh",&_demtsh);
      _trkana->Branch("demtsm",&_demtsm);
    }
// add branches for other tracks
    _trkana->Branch("uem",&_uemti,TrkInfo::leafnames().c_str());
    _trkana->Branch("dmm",&_dmmti,TrkInfo::leafnames().c_str());
// calorimeter information for the downstream electron track
    _trkana->Branch("demc",&_demc,TrkCaloInfo::leafnames().c_str());
// optionally add MC truth branches
    if(_fillmc){
      _trkana->Branch("demmc",&_demmc,TrkInfoMC::leafnames().c_str());
      _trkana->Branch("demmcgen",&_demmcgen,TrkInfoMCStep::leafnames().c_str());
      _trkana->Branch("demmcent",&_demmcent,TrkInfoMCStep::leafnames().c_str());
      if(_diag > 1)_trkana->Branch("demtshmc",&_demtshmc);
    }
  }

  void TrackAnalysis::analyze(const art::Event& event) {
  // reset
    resetBranches();
  // setup KalDiag.
    if(_fillmc)_kdiag.findMCData(event);
// fill basic event information
    fillEventInfo(event);
// find track collections
    // Get handle to tracks collection
    art::Handle<KalRepPtrCollection> demH;
    event.getByLabel(_demtag,demH);
    KalRepPtrCollection const& demC = *demH;
    // Track-cluster matching
    _cdiag.findData(event);
    // find the best track
    const KalRep* demK = findBestTrack(demC);
    // fill the standard diagnostics
    if(demK != 0){
      _kdiag.fillTrkInfo(demK,_demti);
    // fill calorimeter information. First find the best matching cluster
      if(_cdiag.caloMatchHandle().isValid()){
	TrackClusterMatchCollection const& tcmc = *_cdiag.caloMatchHandle();
	TrackClusterMatchCollection::const_iterator itcm = tcmc.end();
	findBestClusterMatch(tcmc,demK,itcm);
	if(itcm != tcmc.end())
	  _cdiag.fillCaloInfo(*itcm,_demc);
      }
    }
    // fill mC info associated with this track
    if(_fillmc) fillMCInfo(demK);
    // fill this row in the TTree
    _trkana->Fill();
  }

  const KalRep* TrackAnalysis::findBestTrack(KalRepPtrCollection const& kcol) {
// if there aren't any tracks return 0
    const KalRep* krep(0);
    if(kcol.size() > 0){
// for now return the 1st track.  This should pick out the most signal-like track FIXME!!
// There should also be information stored about overlaps in this list
      krep = kcol.front().get();
    }
    return krep;
  }

  void TrackAnalysis::fillMCInfo(const KalRep* demK) {
    art::Ptr<SimParticle> demSP;
    if(demK != 0) {
      _kdiag.findMCTrk(demK,demSP);
     } else if(_kdiag.mcData()._simparts != 0) { 
      // assume the 1st primary particle is the CE.  FIXME!!!
      for ( auto isp = _kdiag.mcData()._simparts->begin(); isp != _kdiag.mcData()._simparts->end(); ++isp ){
	if(isp->second.isPrimary()){
	  demSP = art::Ptr<SimParticle>(_kdiag.mcData()._simparthandle,isp->second.id().asInt());
	  break;
	}
      }
     }
    if(demSP.isNonnull()){
      _kdiag.fillTrkInfoMC(demSP,demK,_demmc);
      _kdiag.fillTrkInfoMCStep(demSP,_demmcgen);
      // find where the particle crosses the tracker entranc
      std::vector<MCStepItr> steps;
      _kdiag.findMCSteps(_kdiag.mcData()._mcvdsteps,demSP->id(),_kdiag.VDids(KalDiag::trackerEnt),steps);
      // if there are 2 crossings, take the 2nd
      if(steps.size() >=2 )
	_kdiag.fillTrkInfoMCStep(steps[1],_demmcent);
      else if(steps.size() > 0)
	_kdiag.fillTrkInfoMCStep(steps[0],_demmcent);
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
    for ( const auto& ievtWt : _evtWttags ) {
      _einfo._evtwt *= event.getValidHandle<EventWeight>( ievtWt )->weight();
    }
    // generator weight
    art::Handle<EventWeight> genWtHandle;
    event.getByLabel(_genWttag, genWtHandle);
    if(genWtHandle.isValid())
      _einfo._genwt = genWtHandle->weight();
    // proton bunch weight
    art::Handle<EventWeight> beamWtHandle;
    event.getByLabel(_beamWttag, beamWtHandle);
    if(beamWtHandle.isValid())
      _einfo._beamwt = beamWtHandle->weight();
    // actual number of protons on target
    art::Handle<ProtonBunchIntensity> PBIHandle;
    event.getByLabel(_PBItag, PBIHandle);
    if(PBIHandle.isValid())
      _einfo._nprotons = PBIHandle->intensity();
  }

  void TrackAnalysis::findBestClusterMatch(TrackClusterMatchCollection const& tcmc,
  const KalRep* krep, 
  TrackClusterMatchCollection::const_iterator& itcm) {
  // for now pick highest-energy cluster.  This should match to the track or sum, FIXME!!
    double emin(-1.0);
    for( auto jtcm= tcmc.begin(); jtcm != tcmc.end(); ++jtcm ) {
      if(jtcm->textrapol()->trk().get() == krep){
	if(jtcm->caloCluster()->energyDep() > emin){
	  itcm = jtcm;
	  emin = jtcm->caloCluster()->energyDep(); 
	}
      }
    }
  }

  void TrackAnalysis::resetBranches() {
    _demti.reset();
    _uemti.reset();
    _dmmti.reset();
    _demc.reset();
    _demmc.reset();
    _uemmc.reset();
    _dmmmc.reset();
    _demmcgen.reset();
    _demmcent.reset();
    // clear vectors
    _demtsh.clear();
    _demtsm.clear();
    _demtshmc.clear();
  }
}  // end namespace mu2e

// Part of the magic that makes this class a module.
// create an instance of the module.  It also registers
using mu2e::TrackAnalysis;
DEFINE_ART_MODULE(TrackAnalysis);
