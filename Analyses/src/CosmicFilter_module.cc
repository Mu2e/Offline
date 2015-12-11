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
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
//#include <boost/shared_ptr.hpp>
#include "fhiclcpp/ParameterSet.h"
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "TrackerGeom/inc/Tracker.hh"

#include "BTrk/TrkBase/HelixParams.hh"
#include "RecoDataProducts/inc/KalRepCollection.hh"

using namespace std;

namespace mu2e {

  class CosmicFilter : public art::EDFilter {
  public:

    struct Hist_t {
      TH1F*  fNTracks;
      TH1F*  fNGoodTracks;
      TH1F*  fTrackD0;
      TH1F*  fTrackZ0[2];
    } fHist;
					// module parameters
    std::string  fTrkPatRecModuleLabel;
    int          fDiagLevel;
    double       fMaxD0;
    double       fMaxZ0;
					// otehr variables
    int          fNTracks;
					// Pointers to histograms & ntuples

  public:

    explicit CosmicFilter(fhicl::ParameterSet const& pset):
      fTrkPatRecModuleLabel(pset.get<std::string>   ("trkPatRecModuleLabel","TrkPatRec")),
      fDiagLevel           (pset.get<int>           ("diagLevel",   0 )),
      fMaxD0               (pset.get<double>        ("maxD0"    , 200.)), // in mm
      fMaxZ0               (pset.get<double>        ("maxZ0"    ,1000.))  // in mm
    {
    }

    virtual ~CosmicFilter() {}

    virtual void beginJob();
    virtual bool filter  (art::Event& event);

  };

//-----------------------------------------------------------------------------
  void CosmicFilter::beginJob(){

    // Get access to the TFile service.
    art::ServiceHandle<art::TFileService> tfs;

    // book histograms.
    fHist.fTrackD0     = tfs->make<TH1F>("trk_d0"  ,"Track D0"      ,500, -500., 500.);
    fHist.fNTracks     = tfs->make<TH1F>("ntrk"    ,"N(tracks)"     ,100,    0., 100.);
    fHist.fNGoodTracks = tfs->make<TH1F>("ngtrk"   ,"N(good tracks)",100,    0., 100.);
    fHist.fTrackZ0[0]  = tfs->make<TH1F>("trk_z0_0","Track Z0[0]"   ,500,-2500.,2500.);
    fHist.fTrackZ0[1]  = tfs->make<TH1F>("trk_z0_1","Track Z0[1]"   ,500,-2500.,2500.);
  }

//-----------------------------------------------------------------------------
  bool CosmicFilter::filter(art::Event& anEvent) {

    const KalRep*   trk;
    double          d0, z0;
    int             n_good_tracks;

    bool rc = false;

    art::Handle<KalRepCollection> krepsHandle;
    anEvent.getByLabel(fTrkPatRecModuleLabel,"DownstreameMinus", krepsHandle);
    const KalRepCollection*  kreps = krepsHandle.product();

    fNTracks = kreps->size();

    n_good_tracks = 0;
    for (int i=0; i<fNTracks; i++) {
      trk = kreps->get(i);

      d0 = trk->helix(0).d0();
      z0 = trk->helix(0).z0();

      fHist.fTrackD0->Fill(d0);
      fHist.fTrackZ0[0]->Fill(z0);
      if (fabs(d0) < fMaxD0) fHist.fTrackZ0[1]->Fill(z0);

      if ((fabs(d0) < fMaxD0) && (fabs(z0) < fMaxZ0)) {
	n_good_tracks += 1;
      }
    }

    fHist.fNTracks->Fill(fNTracks);
    fHist.fNGoodTracks->Fill(n_good_tracks);

    rc = (n_good_tracks > 0);

    return rc;
  }
}

using mu2e::CosmicFilter;
DEFINE_ART_MODULE(CosmicFilter);
