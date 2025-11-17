///////////////////////////////////////////////////////////////////////////////
// filter out cosmics events
// -------------------------
// 2013-03-01: first day of the sequester - require an event to have:
//
// 1. a reconstructed track with |d0| < 150mm
//
///////////////////////////////////////////////////////////////////////////////

// C++ includes
#include <iostream>
#include <stdexcept>
#include <string>

#include "TH1.h"

// Framework includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
//#include <boost/shared_ptr.hpp>
#include "fhiclcpp/ParameterSet.h"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

#include "BTrkLegacy/inc/HelixParams.hh"
#include "Offline/RecoDataProducts/inc/KalRepPtrCollection.hh"

using namespace std;

namespace mu2e {

  class CosmicFilter : public art::EDFilter {
  public:

    struct Hist_t {
      TH1F*  fNTracks;                        // number of reconstructed tracks
      TH1F*  fNGoodTracks;                // N(tracks passing all cuts)
      TH1F*  fTrackD0;                        // impact parameter
      TH1F*  fTrackZ0[2];                // Z0 in the point of closest approach
      TH1F*  fP;                         // track momentum
    } _hist;
                                        // module parameters
    std::string  fTrkPatRecModuleLabel;
    int          fDiagLevel;
    double       fMaxD0;
    double       fMaxZ0;
                                        // otehr variables
    int          _nTracks;
                                        // Pointers to histograms & ntuples

  public:

    explicit CosmicFilter(fhicl::ParameterSet const& pset);
    virtual ~CosmicFilter() {}

    virtual void beginJob();
    virtual bool filter  (art::Event& event);

  };

    CosmicFilter::CosmicFilter(fhicl::ParameterSet const& pset):
    EDFilter{pset},
    fTrkPatRecModuleLabel(pset.get<std::string>   ("trkPatRecModuleLabel")),
    fDiagLevel           (pset.get<int>           ("diagLevel"           )),
    fMaxD0               (pset.get<double>        ("maxD0"               )), // in mm
    fMaxZ0               (pset.get<double>        ("maxZ0"               ))  // in mm
    {}

//-----------------------------------------------------------------------------
// Get access to the TFile service and book histograms
//-----------------------------------------------------------------------------
  void CosmicFilter::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;

    _hist.fTrackD0     = tfs->make<TH1F>("trk_d0"  ,"Track D0"       ,500, -500., 500.);
    _hist.fNTracks     = tfs->make<TH1F>("ntrk"    ,"N(tracks)"      ,100,    0., 100.);
    _hist.fNGoodTracks = tfs->make<TH1F>("ngtrk"   ,"N(good tracks)" ,100,    0., 100.);
    _hist.fTrackZ0[0]  = tfs->make<TH1F>("trk_z0_0","Track Z0[0]"    ,500,-2500.,2500.);
    _hist.fTrackZ0[1]  = tfs->make<TH1F>("trk_z0_1","Track Z0[1]"    ,500,-2500.,2500.);
    _hist.fP           = tfs->make<TH1F>("p"       ,"Track Mom (S=0)",200,    0., 400.);
  }

//-----------------------------------------------------------------------------
  bool CosmicFilter::filter(art::Event& anEvent) {

    const KalRep*   trk;
    double          d0, z0;
    int             n_good_tracks;

    bool rc = false;

    art::Handle<KalRepPtrCollection> krepsHandle;
    anEvent.getByLabel(fTrkPatRecModuleLabel,"DownstreameMinus", krepsHandle);
    const KalRepPtrCollection*  list_of_kreps(0);


    _nTracks = 0;
    if (krepsHandle.isValid()) {
      list_of_kreps = krepsHandle.product();
      _nTracks      = list_of_kreps->size();
    }


    n_good_tracks = 0;
    for (int i=0; i<_nTracks; i++) {
      trk = list_of_kreps->at(i).get();

      d0 = trk->helix(0).d0();
      z0 = trk->helix(0).z0();

      _hist.fTrackD0->Fill(d0);
      _hist.fTrackZ0[0]->Fill(z0);
      if (fabs(d0) < fMaxD0) _hist.fTrackZ0[1]->Fill(z0);

      if ((fabs(d0) < fMaxD0) && (fabs(z0) < fMaxZ0)) {
        n_good_tracks += 1;
      }

      CLHEP::Hep3Vector mom = trk->momentum(0); // at S=0
      _hist.fP->Fill(mom.mag());
    }

    _hist.fNTracks->Fill(_nTracks);
    _hist.fNGoodTracks->Fill(n_good_tracks);

    rc = (n_good_tracks > 0);

    return rc;
  }
}

using mu2e::CosmicFilter;
DEFINE_ART_MODULE(CosmicFilter)
