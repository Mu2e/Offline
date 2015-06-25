//
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
#include "MCDataProducts/inc/EventWeight.hh"
#include "Mu2eUtilities/inc/SimpleSpectrum.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"

// ROOT incldues
#include "TH1F.h"

// Need this for the BaBar headers.
using namespace CLHEP;

// BaBar includes
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
// mu2e tracking
#include "KalmanTests/inc/TrkFitDirection.hh"
#include "KalmanTests/inc/KalDiag.hh"

// C++ includes.
#include <iostream>
#include <string>

//G4Beamline includes
#include "MCDataProducts/inc/G4BeamlineInfoCollection.hh"

// This is fragile and needs to be last until CLHEP is
// properly qualified and included in the BaBar classes.
#include "KalmanTests/inc/KalRepCollection.hh"

using namespace std;

namespace mu2e {

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
    vector<art::InputTag> _evtWtModules;
    TrkParticle _tpart;
    TrkFitDirection _fdir;
    std::string _iname;

    bool haveG4BL;// = g4beamlineData.isValid();
    // diagnostic of Kalman fit
    KalDiag _kdiag;

    // Control level of printout.
    int _verbosity;
    int _maxPrint;
    // whether or not to include MC info for empty events
    bool _processEmpty;

    // Histograms
    TH1F* _hNTracks;
    TH1F* _hfitCL;
    TH1F* _hChisq;
    TH1F* _hmomentum0;
    TH1F* _hdp;
    TH1F* _hz0;

    TTree* _trkdiag;


    Int_t _trkid,_eventid, _runid, _subrunid;
    Double_t _evtwt; 
    Float_t g4bl_weight;

  };

  ReadKalFits::ReadKalFits(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset),
    _fitterModuleLabel(pset.get<string>("fitterModuleLabel")),
    _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
    _evtWtModules( pset.get<std::vector<art::InputTag>>("eventWeightModules",std::vector<art::InputTag>() ) ),
    _tpart((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
    _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
    _kdiag(pset.get<fhicl::ParameterSet>("KalDiag",fhicl::ParameterSet())),
    _verbosity(pset.get<int>("verbosity",0)),
    _maxPrint(pset.get<int>("maxPrint",0)),
    _processEmpty(pset.get<bool>("processEmpty",true)),
    _hNTracks(0),
    _hfitCL(0),
    _hChisq(0),
    _hmomentum0(0),
    _hdp(0),
    _hz0(0),
    _trkdiag(0) {
// construct the data product instance name
    _iname = _fdir.name() + _tpart.name();
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
    _trkdiag->Branch("evtwt",&_evtwt,"evtwt/d");
    _trkdiag->Branch("g4bl_weight",&g4bl_weight,"g4bl_weight/f");
}

  // For each event, look at tracker hits and calorimeter hits.
  void ReadKalFits::analyze(const art::Event& event) {
    //    cout << "Enter ReadKalFits:: analyze: " << _verbosity << endl;

    _eventid = event.event();
    _runid = event.run();
    _subrunid = event.subRun();
    if(!_kdiag.findMCData(event)){
      throw cet::exception("RECO")<<"mu2e::KalDiag: MC data missing or incoplete" << std::endl;
    }
    // Get handle to calorimeter hit collection.
    art::Handle<KalRepCollection> trksHandle;
    event.getByLabel(_fitterModuleLabel,_iname,trksHandle);
    KalRepCollection const& trks = *trksHandle;

    art::Handle<G4BeamlineInfoCollection> g4beamlineData;
    event.getByLabel(_generatorModuleLabel, g4beamlineData);
   
    if ( _verbosity > 0 && _eventid <= _maxPrint ){
      cout << "ReadKalmanFits  for event: " << event.id() << "  Number of fitted tracks: " << trks.size() << endl;
    }

    // Modify event weight
    _evtwt = 1.;
    for ( const auto& ievtWt : _evtWtModules ) {
      _evtwt *= event.getValidHandle<EventWeight>( ievtWt )->weight();
    }
  //	g4bl_weight=1;
    haveG4BL = g4beamlineData.isValid();
    if ( haveG4BL ) haveG4BL = (g4beamlineData->size()==1);
     if( haveG4BL ) { 
      G4BeamlineInfo const& extra = g4beamlineData->at(0);
      g4bl_weight=extra.weight();
    } else{
      g4bl_weight=1;
    }

     _hNTracks->Fill( trks.size() );
     _trkid = -1;

    for ( size_t i=0; i< trks.size(); ++i ){
      _trkid = i;
      KalRep const* krep = trks.get(i);
      if ( !krep ) continue;

     _kdiag.kalDiag(krep);

      // For some quantities you require the concrete representation, not
      // just the base class.
      // Fill a histogram.
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
             << i            << " "
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
    // if there are no tracks, enter dummies
    if(trks.size() == 0 && _processEmpty){
      _kdiag.kalDiag(0);
    }
  }

}  // end namespace mu2e

// Part of the magic that makes this class a module.
// create an instance of the module.  It also registers
using mu2e::ReadKalFits;
DEFINE_ART_MODULE(ReadKalFits);
