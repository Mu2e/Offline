//
// Read the tracks added to the event by KalFitTest_module.
//
// $Id: ReadKalFits_module.cc,v 1.9 2012/07/20 22:38:12 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/20 22:38:12 $
//
// Original author Rob Kutschke
//

// Mu2e includes.

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
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkRep.hh"
#include "KalmanTrack/KalRep.hh"
#include "KalmanTests/inc/KalFitMC.hh"

// C++ includes.
#include <iostream>
#include <string>

// This is fragile and needs to be last until CLHEP is
// properly qualified and included in the BaBar classes.
#include "KalmanTests/inc/TrkRecoTrkCollection.hh"

using namespace std;

namespace mu2e {

  class ReadKalFits : public art::EDAnalyzer {
  public:

    explicit ReadKalFits(fhicl::ParameterSet const& pset);
    virtual ~ReadKalFits() { }

    void beginJob();
    void analyze(const art::Event& e);

    // DIO spectrum
    static double DIOspectrum(double eenergy);

  private:

    // Module label of the module that performed the fits.
    std::string _fitterModuleLabel;
    // whether to weight the DIO or not
    bool _weight;
    // diagnostic of Kalman fit
    KalFitMC _kfitmc;

    // Control level of printout.
    int _verbosity;
    int _maxPrint;

    // Histograms
    TH1F* _hNTracks;
    TH1F* _hfitCL;
    TH1F* _hChisq;
    TH1F* _hmomentum0;
    TH1F* _hdp;
    TH1F* _hz0;

    TTree* _trkdiag;


    Int_t _trkid,_eventid;
    Float_t _diowt;

  };

  ReadKalFits::ReadKalFits(fhicl::ParameterSet const& pset):
    _fitterModuleLabel(pset.get<string>("fitterModuleLabel")),
    _weight(pset.get<bool>("WeightEvents",true)),
    _kfitmc(pset.get<fhicl::ParameterSet>("KalFitMC")),
    _verbosity(pset.get<int>("verbosity",0)),
    _maxPrint(pset.get<int>("maxPrint",0)),
    _hNTracks(0),
    _hfitCL(0),
    _hChisq(0),
    _hmomentum0(0),
    _hdp(0),
    _hz0(0),
    _trkdiag(0){
  }

  void ReadKalFits::beginJob( ){
    art::ServiceHandle<art::TFileService> tfs;
    _hNTracks   = tfs->make<TH1F>( "hNTracks",  "Number of tracks per event.",         10,     0.,   10. );
    _hfitCL     = tfs->make<TH1F>( "hfitCL",    "Confidence Level of the Track fit.",  50,     0.,    1. );
    _hChisq     = tfs->make<TH1F>( "hChisq",    "Chisquared of the Track fit.",       100,     0.,  500. );
    _hmomentum0 = tfs->make<TH1F>( "hmomentum", "Reco Momentum at front face",        100,    70.,  110. );
    _hdp        = tfs->make<TH1F>( "hdp",       "Momentum Change across Tracker",     100,   -10.,    0. );
    _hz0        = tfs->make<TH1F>( "hz0",       "z of start valid region (cm)",       100, -1500., -800. );
    _trkdiag    = _kfitmc.createTrkDiag();
// add local branches
    _trkdiag->Branch("eventid",&_eventid,"eventid/I");
    _trkdiag->Branch("trkid",&_trkid,"trkid/I");
    _trkdiag->Branch("diowt",&_diowt,"diowt/f");
    _eventid = 0;
  }

  // For each event, look at tracker hits and calorimeter hits.
  void ReadKalFits::analyze(const art::Event& event) {
    cout << "Enter ReadKalFits:: analyze: " << _verbosity << endl;

    _eventid++;
    _kfitmc.findMCData(event);
    // DIO spectrum weight; use the generated momenum magnitude
    double ee = _kfitmc._mcmom;
    _diowt = 1.0;
    if(_weight)
      _diowt = DIOspectrum(ee);
    // Get handle to calorimeter hit collection.
    art::Handle<TrkRecoTrkCollection> trksHandle;
    event.getByLabel(_fitterModuleLabel,trksHandle);
    TrkRecoTrkCollection const& trks = *trksHandle;

    if ( _verbosity > 0 && _eventid <= _maxPrint ){
      cout << "ReadKalmanFits  for event: " << event.id() << "  Number of fitted tracks: " << trks.size() << endl;
    }

    _hNTracks->Fill( trks.size() );
    _trkid = -1;
    for ( size_t i=0; i< trks.size(); ++i ){
      _trkid = i;
      TrkRecoTrk const& trk = *trks[i];
      TrkRep const* trep = trk.getRep(PdtPid::electron);
      if ( !trep ) continue;

      _kfitmc.trkDiag(trk);

      // For some quantities you require the concrete representation, not
      // just the base class.
      KalRep const* krep = dynamic_cast<KalRep const*>(trep);
      if ( !krep ) continue;

      // Fill a histogram.
      _hfitCL->Fill(krep->chisqConsistency().likelihood() );
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
    if(trks.size() == 0){
      TrkKalFit dummy;
      _kfitmc.trkDiag(dummy);
    }
  }


  double ReadKalFits::DIOspectrum(double ee) {
    double mal(25133);
//    double mmu(105.654);
    double emu(105.194);
//    double emue(104.973);
//    double me(0.511);
    double a5(8.6434e-17);
    double a6(1.16874e-17);
    double a7(-1.87828e-19);
    double a8(9.16327e-20);
    double delta = emu - ee - ee*ee/(2*mal);
    return a5*pow(delta,5) + a6*pow(delta,6) + a7*pow(delta,7) + a8*pow(delta,8);
  }

}  // end namespace mu2e

// Part of the magic that makes this class a module.
// create an instance of the module.  It also registers
using mu2e::ReadKalFits;
DEFINE_ART_MODULE(ReadKalFits);
