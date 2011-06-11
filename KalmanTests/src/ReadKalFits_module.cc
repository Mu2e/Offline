//
// Read the tracks added to the event by KalFitTest_module.
//
// $Id: ReadKalFits_module.cc,v 1.1 2011/06/11 03:19:08 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/11 03:19:08 $
//
// Original author Rob Kutschke
//

// Mu2e includes.

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
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

  private:

    // Module label of the module that performed the fits.
    std::string _fitterModuleLabel;

    // Histograms
    TH1F* _hNTracks;
    TH1F* _hfitCL;
    TH1F* _hChisq;

  };

  ReadKalFits::ReadKalFits(fhicl::ParameterSet const& pset):
    _fitterModuleLabel(pset.get<string>("fitterModuleLabel")),
    _hNTracks(0),
    _hfitCL(0),
    _hChisq(0){
  }

  void ReadKalFits::beginJob( ){
    art::ServiceHandle<art::TFileService> tfs;
    _hNTracks = tfs->make<TH1F>( "hNTracks", "Number of tracks per event.",         10, 0.,   10. );
    _hfitCL   = tfs->make<TH1F>( "hfitCL",   "Confidence Level of the Track fit.",  50, 0.,    1. );
    _hChisq   = tfs->make<TH1F>( "hChisq",   "Chisquared of the Track fit.",       100, 0.,  500. );
  }

  // For each event, look at tracker hits and calorimeter hits.
  void ReadKalFits::analyze(const art::Event& event) {

    // Get handle to calorimeter hit collection.
    art::Handle<TrkRecoTrkCollection> trksHandle;
    event.getByLabel(_fitterModuleLabel,trksHandle);
    TrkRecoTrkCollection const& trks = *trksHandle;

    _hNTracks->Fill( trks.size() );

    for ( size_t i=0; i< trks.size(); ++i ){

      TrkRecoTrk const& trk = trks[i];
      TrkRep const* trep = trk.getRep(PdtPid::electron);
      if ( !trep ) continue;

      // For some quantities you require the concrete representation, not
      // just the base class.
      KalRep const* krep = dynamic_cast<KalRep const*>(trep);
      if ( !krep ) continue;

      // Fill a histogram.
      _hfitCL->Fill(trep->chisqConsistency().likelihood() );
      _hChisq->Fill(trep->chisqConsistency().chisqValue() );

    }
  }

}  // end namespace mu2e

// Part of the magic that makes this class a module.
// create an instance of the module.  It also registers
using mu2e::ReadKalFits;
DEFINE_ART_MODULE(ReadKalFits);
