// Read the tracks added to the event by KalFitTest_module.
//
// $Id: ReadKalFits_module.cc,v 1.29 2014/08/28 19:26:04 brownd Exp $
// $Author: brownd $
// $Date: 2014/08/28 19:26:04 $
//
// Original author Rob Kutschke
//

// Mu2e includes
#include "GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "MCDataProducts/inc/EventWeight.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"

// ROOT incldues
#include "TH1F.h"
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

  class ReadKalFits : public art::EDAnalyzer {

  public:

    explicit ReadKalFits(fhicl::ParameterSet const& pset);
    virtual ~ReadKalFits() { }

    void beginJob();
    void analyze(const art::Event& e);

  private:

    typedef vector<string> VS;

    // Module label of the module that performed the fits.
    std::string _fitterModuleLabel;
    // Label of the generator.
    std::string _generatorModuleLabel;
    // Label of the event-weighting module
    art::InputTag _genWtModule;
    art::InputTag _beamWtModule;
    art::InputTag _PBIModule;
    vector<art::InputTag> _evtWtModules;
    TrkParticle _tpart;
    TrkFitDirection _fdir;
    // diagnostic of Kalman fit
    KalDiag _kdiag;
    // calorimeter diagnostics.  Unfortunately the PID objects have no default constructor
    // so this MUST be created on the heap to maintain backwards compatibility.  FIXME!!
    TrkCaloDiag* _cdiag;
    // Control level of printout.
    int _verbosity;
    int _maxPrint;
    // whether or not to include MC info for empty events
    bool _processEmpty;
    // whether or not to include calorimeter information
    bool _addCalo;
// event data
    art::Handle<KalRepPtrCollection> _trksHandle;
    art::Handle<TrackClusterMatchCollection> _caloMatchHandle;

    // Histograms
    TH1F* _hNTracks;
    TH1F* _hfitCL;
    TH1F* _hChisq;
    TH1F* _hmomentum0;
    TH1F* _hdp;
    TH1F* _hz0;
// main diagnostic TTree, based on KalDiag
    TTree* _trkdiag;
//  local branches
    Int_t _trkid,_eventid, _runid, _subrunid;
    Double_t _evtwt, _beamwt, _genwt;
    Int_t _nprotons;
    Float_t _g4bl_weight;
    Int_t _ntrks, _ntshared;
    std::vector<TrkHitShare> _overlaps;

    // helper functions
    void countSharedHits(KalRepPtrCollection const& trks,
	std::vector<TrkHitShare>& overlaps);
    void findWeights(const art::Event& event);
  };

  ReadKalFits::ReadKalFits(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset),
    _fitterModuleLabel(pset.get<string>("fitterModuleLabel")),
    _generatorModuleLabel(pset.get<std::string>("generatorModule", "generate")),
    _genWtModule( pset.get<art::InputTag>("generatorWeightModule",art::InputTag()) ),
    _beamWtModule( pset.get<art::InputTag>("beamWeightModule",art::InputTag()) ),
    _PBIModule( pset.get<art::InputTag>("ProtonBunchIntensityModule",art::InputTag("protonBunchIntensity")) ),
    _evtWtModules( pset.get<std::vector<art::InputTag>>("eventWeightModules",std::vector<art::InputTag>() ) ),
    _tpart((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
    _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
    _kdiag(pset.get<fhicl::ParameterSet>("KalDiag",fhicl::ParameterSet())),
    _cdiag(0),
    _verbosity(pset.get<int>("verbosity",0)),
    _maxPrint(pset.get<int>("maxPrint",0)),
    _processEmpty(pset.get<bool>("processEmpty",true)),
    _addCalo(pset.get<bool>("addCalo",false)),
    _hNTracks(0),
    _hfitCL(0),
    _hChisq(0),
    _hmomentum0(0),
    _hdp(0),
    _hz0(0),
    _trkdiag(0)
  {
// construct the data product instance name
    if(_addCalo){
      _cdiag = new TrkCaloDiag(_tpart,_fdir,pset.get<fhicl::ParameterSet>("TrkCaloDiag",fhicl::ParameterSet()));
    }
  }

  void ReadKalFits::beginJob( ){
    art::ServiceHandle<art::TFileService> tfs;
    _hNTracks   = tfs->make<TH1F>( "hNTracks",  "Number of tracks per event.",         10,     0.,   10. );
    _hfitCL     = tfs->make<TH1F>( "hfitCL",    "Confidence Level of the Track fit.",  50,     0.,    1. );
    _hChisq     = tfs->make<TH1F>( "hChisq",    "Chisquared of the Track fit.",       100,     0.,  500. );
    _hmomentum0 = tfs->make<TH1F>( "hmomentum", "Reco Momentum at front face",        100,    70.,  110. );
    _hdp        = tfs->make<TH1F>( "hdp",       "Momentum Change across Tracker",     100,   -10.,    0. );
    _hz0        = tfs->make<TH1F>( "hz0",       "z of start valid region (cm)",       100, -1500., -800. );
    _trkdiag    = _kdiag.createTrkDiag();
    // add local branches
    _trkdiag->Branch("eventid",&_eventid,"eventid/I");
    _trkdiag->Branch("runid",&_runid,"runid/I");
    _trkdiag->Branch("subrunid",&_subrunid,"subrunid/I");
    _trkdiag->Branch("trkid",&_trkid,"trkid/I");
    _trkdiag->Branch("genwt",&_genwt,"genwt/d");
    _trkdiag->Branch("beamwt",&_beamwt,"beamwt/d");
    _trkdiag->Branch("evtwt",&_evtwt,"evtwt/d");
    _trkdiag->Branch("nprotons",&_nprotons,"nprotons/I");
    _trkdiag->Branch("g4bl_weight",&_g4bl_weight,"g4bl_weight/f");
    _trkdiag->Branch("ntrks",&_ntrks,"ntrks/I");
    _trkdiag->Branch("ntshared",&_ntshared,"ntshared/I");
    _trkdiag->Branch("overlaps",&_overlaps);
    if(_addCalo){
      _cdiag->addBranches(_trkdiag);
    }
  }

  // For each event, look at tracker hits and calorimeter hits.
  void ReadKalFits::analyze(const art::Event& event) {
// basic event information
    _eventid = event.event();
    _runid = event.run();
    _subrunid = event.subRun();
// fill event weight information
    findWeights(event);
// collect all the event information
    if(!_kdiag.findMCData(event))
      throw cet::exception("RECO")<<"mu2e::ReadKalFits: MC data missing or incomplete" << std::endl;
    // Get handle to tracks collection
    art::Handle<KalRepPtrCollection> trksHandle;
    event.getByLabel(_fitterModuleLabel,trksHandle);
    KalRepPtrCollection const& trks = *trksHandle;
    // Track-cluster matching
    if(_addCalo)_cdiag->findData(event);
    _hNTracks->Fill( trks.size() );
    // initialize counting variables
     _ntrks = trks.size();
     _ntshared = -1;
     _trkid = -1;
    // search for tracks which share hits
    _overlaps.clear();
    countSharedHits(trks,_overlaps);
    // diagnostic printout
    if ( _verbosity > 0 && _eventid <= _maxPrint ){
      cout << "ReadKalmanFits  for event: " << event.id() << "  Number of fitted tracks: " << trks.size() << endl;
    }
    // if there are no tracks, make an entry for a 'null' track.  This keeps the MC bookkeeping complete
    if(trks.size() == 0 && _processEmpty){
      _kdiag.kalDiag(0);
      if(_addCalo)_cdiag->addCaloInfo(0);
    }
// main loop over tracks
    for ( size_t itrk=0; itrk< trks.size(); ++itrk ){
    // we don't have a real TrackID in Mu2e: just use the list index
      _trkid = itrk;
      KalRep const* krep = trks.at(itrk).get();
      if ( !krep ) continue;
      _ntshared = 0;
      for(const TrkHitShare& ihs: _overlaps) {
	if(ihs._trk1 == itrk || ihs._trk2 == itrk)
	  ++_ntshared;
      }
      // if requested, find matching calorimeter information; can be more than 1
      if(_addCalo)_cdiag->addCaloInfo(krep);
      // fill the standard diagnostics
      _kdiag.kalDiag(krep);


      // Fill some simple histograms
      _hfitCL->Fill(krep->chisqConsistency().significanceLevel() );
      _hChisq->Fill(krep->chisqConsistency().chisqValue() );

      double s0   = krep->startValidRange();
      double sEnd = krep->endValidRange();

      Hep3Vector momentum0   = krep->momentum(s0);
      Hep3Vector momentumEnd = krep->momentum(sEnd);
      HepPoint   pos0        = krep->position(s0);
      double dp              = momentumEnd.mag()-momentum0.mag();

      _hmomentum0->Fill(momentum0.mag());
      _hdp->Fill(dp);

      _hz0->Fill(pos0.z());

      if ( _verbosity > 1 && _eventid <= _maxPrint ){
        cout << "   Fitted track: "
             << itrk            << " "
             << krep->nDof() << " "
             << krep->chisqConsistency().likelihood()  << " | "
             << s0                        << " "
             << krep->startFoundRange()   << " | "
             << sEnd                      << " "
             << krep->endFoundRange()     << " | "
             << momentum0.mag()           << " "
             << momentumEnd.mag()         << " | "
             << pos0.z()
             << endl;
      }

    }
  }

  void ReadKalFits::countSharedHits(KalRepPtrCollection const& trks,
    std::vector<TrkHitShare>& overlaps) {
    if(trks.size() > 1){
      for(size_t itrk =0; itrk < trks.size(); ++itrk){
	KalRep const* ikrep = trks.at(itrk).get();
	TrkStrawHitVector ihits;
	convert(ikrep->hitVector(),ihits);
	for(size_t jtrk = itrk++; jtrk < trks.size(); ++jtrk ) {
	  KalRep const* jkrep = trks.at(jtrk).get();
	  TrkStrawHitVector jhits;
	  convert(jkrep->hitVector(),jhits);
	  unsigned nhshared(0);
	  for(const TrkStrawHit* ihit: ihits){
	    for(const TrkStrawHit* jhit: jhits){
	      if(ihit->isActive() && jhit->isActive() && &ihit->comboHit() == &jhit->comboHit())
		++nhshared;
	    }
	  }
	  if(nhshared > 0){
	    TrkHitShare share;
	    share._nhshared = nhshared;
	    if(ikrep->nActive() > jkrep->nActive()){
	      share._trk1 = itrk;
	      share._trk2 = jtrk;
	      share._f1 = nhshared/float(ikrep->nActive());
	      share._f2 = nhshared/float(jkrep->nActive());
	    } else {
	      share._trk2 = itrk;
	      share._trk1 = jtrk;
	      share._f1 = nhshared/float(jkrep->nActive());
	      share._f2 = nhshared/float(ikrep->nActive());
	    }
	    overlaps.push_back(share);
	  }
	}
      }
    }
  }

  void ReadKalFits::findWeights( const art::Event& event) {
   
    // get event weight product
    _genwt = _beamwt = _evtwt = _g4bl_weight = 1.; 
    _nprotons=-1;
    // total weight is the product of all weights
    for ( const auto& ievtWt : _evtWtModules ) {
      _evtwt *= event.getValidHandle<EventWeight>( ievtWt )->weight();
    }
    // generator weight
    art::Handle<EventWeight> genWtHandle;
    event.getByLabel(_genWtModule, genWtHandle);
    if(genWtHandle.isValid())
      _genwt = genWtHandle->weight();
    // proton bunch weight
    art::Handle<EventWeight> beamWtHandle;
    event.getByLabel(_beamWtModule, beamWtHandle);
    if(beamWtHandle.isValid())
      _beamwt = beamWtHandle->weight();
    // actual number of protons on target
    art::Handle<ProtonBunchIntensity> PBIHandle;
    event.getByLabel(_PBIModule, PBIHandle);
    if(PBIHandle.isValid())
      _nprotons = PBIHandle->intensity();
    // g4beamline
    art::Handle<G4BeamlineInfoCollection> g4beamlineData;
    event.getByLabel(_generatorModuleLabel, g4beamlineData);
    if( g4beamlineData.isValid() && g4beamlineData->size()==1) {
      G4BeamlineInfo const& extra = g4beamlineData->at(0);
      _g4bl_weight=extra.weight();
    }
  }


}  // end namespace mu2e

// Part of the magic that makes this class a module.
// create an instance of the module.  It also registers
using mu2e::ReadKalFits;
DEFINE_ART_MODULE(ReadKalFits);
