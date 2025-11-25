// An art analyzer module for the initial calibration of interdetector timing.

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "TTree.h"
#include "art_root_io/TFileService.h"
#include <iostream>
#include <vector>

using namespace std;

namespace mu2e {

  class InterdetectorTimeBootstrap : public art::EDAnalyzer {

     public:

      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<int> diagLevel{Name("diagLevel"),Comment("diag level"),0};
        fhicl::Atom<art::InputTag> kalTag{Name("kalseedModuleLabel"), Comment("kal seed info")};
        fhicl::Atom<art::InputTag> caloTag{Name("caloClusterModuleLabel"), Comment("calo cluster info")};
        fhicl::Atom<art::InputTag> crvTag{Name("crvClusterModuleLabel"), Comment("crv info")};
        fhicl::Atom<double> maxDeltaR{Name("maxDeltaR"),Comment("maxDeltaR"),1000.};
        fhicl::Atom<double> maxDeltaT{Name("maxDeltaT"),Comment("maxDeltaT"),1000.};
      };
      typedef art::EDAnalyzer::Table<Config> Parameters;

      explicit InterdetectorTimeBootstrap(const Parameters& conf);
      virtual ~InterdetectorTimeBootstrap() {}


      virtual void beginJob();
      virtual void endJob();
      virtual void analyze(const art::Event& e) override;

     private:
      Config _conf;
      int _diagLevel;
      art::InputTag _kalseedToken;
      art::InputTag _caloClusterToken;
      art::InputTag _crvClusterToken;
      double _maxDeltaR;
      double _maxDeltaT;
      TTree* _timingTree;
      double _deltaT_trkcal;
      double _deltaT_trkcrv;
};


  InterdetectorTimeBootstrap::InterdetectorTimeBootstrap(const Parameters& conf):
  art::EDAnalyzer(conf),
  _diagLevel(conf().diagLevel()),
  _kalseedToken(conf().kalTag()),
  _caloClusterToken(conf().caloTag()),
  _crvClusterToken(conf().crvTag()),
  _maxDeltaR(conf().maxDeltaR()),
  _maxDeltaT(conf().maxDeltaT())
  {}

  //InterdetectorTimeBootstrap::~InterdetectorTimeBootstrap() {}

  void InterdetectorTimeBootstrap::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    _timingTree = tfs->make<TTree>("timingTree", "Timing Differences");
    _timingTree->Branch("deltaT_trkcal", &_deltaT_trkcal, "deltaT_trkcal/F");
    _timingTree->Branch("deltaT_trkcrv", &_deltaT_trkcrv, "deltaT_trkcrv/F");
  }

  void InterdetectorTimeBootstrap::analyze(const art::Event& event) {
    // Access event data products for all three main detectors:
    art::Handle<KalSeedCollection> tracksHandle;
    event.getByLabel(_kalseedToken, tracksHandle);
    const KalSeedCollection& tracks = *tracksHandle;

    art::Handle<CaloClusterCollection> clustersHandle;
    event.getByLabel(_caloClusterToken, clustersHandle);
    const CaloClusterCollection& clusters = *clustersHandle;

    art::Handle<CrvCoincidenceClusterCollection> crvHandle;
    event.getByLabel(_crvClusterToken, crvHandle);
    const CrvCoincidenceClusterCollection& CRVclusters = *crvHandle;

    // Get geometry handles
    //GeomHandle<Calorimeter> calo;
    //GeomHandle<Tracker> tracker;
    //GeomHandle<CosmicRayShield> crv;

    for (const auto& track : tracks) {
      // Use the track' parameters for extrapolation
      CLHEP::Hep3Vector track_pos_at_calo; // Position of track extrapolated to the calo face

      // Use specific Mu2e geometry functions to find the intersection
      bool intersects = true; // Use geometry to check intersection

      if (intersects) {
        for (const auto& cluster : clusters) {
          const CLHEP::Hep3Vector& cluster_pos = cluster.cog3Vector();
          // Simple geometric match: check proximity
          if ((track_pos_at_calo - cluster_pos).mag() < _maxDeltaR) {
            // Calculate time difference
            // Time of flight correction would be applied here in a real scenario
            double time_calo = cluster.time();
            // The track time needs to be adjusted for the time of flight from the track origin to the calo
            double time_track_extrapolated = track.time() + (track_pos_at_calo.mag() / 299.792458); // TOF in ns
            _deltaT_trkcal = time_calo - time_track_extrapolated;

            // Fill the NTuple for later analysis (e.g., histogramming the mean)
            _timingTree->Fill();
          }
        }
      }
    }
  }
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::InterdetectorTimeBootstrap)
