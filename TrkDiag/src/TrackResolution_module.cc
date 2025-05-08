//
// Resolve alternate track fit hypotheses
//
// Original author M. MacKenzie (2025)
//

// framework
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art_root_io/TFileService.h"

// data
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"

// ROOT
#include "TFile.h"
#include "TH1.h"

// C++
#include <iostream>
#include <map>

namespace mu2e
{

  class TrackResolution : public art::EDProducer
  {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<art::InputTag> seedTag{Name("TrackClusterCollection"), Comment("Input clustering of tracks")};
      fhicl::Atom<bool> useCaloHit{Name("prioritizeCaloHit"), Comment("Prioritize tracks that include a calorimeter hit"), false};
      fhicl::Atom<bool> makeHists{Name("makeHistograms"), Comment("Make debug histograms"), false};
      fhicl::Atom<int> debug{Name("debugLevel"), Comment("Debug printout level"), 0};
    };

    using Parameters = art::EDProducer::Table<Config>;
    TrackResolution(const Parameters& conf);

  private:
    void produce(art::Event& event) override;
    int resolve(const KalSeed* k_1, const KalSeed* k_2);

    art::InputTag _seedTag;
    bool _useCaloHit;
    bool _makeHists;
    int _debug;

    enum {kFirst, kSecond, kNeither}; // possible resolution results

    struct Hist_t {
      TH1* fitcon;
    };
    Hist_t _hists;
  };

  //--------------------------------------------------------------------------------------------------------
  TrackResolution::TrackResolution(const Parameters& conf) :
    art::EDProducer{conf},
    _seedTag(conf().seedTag()),
    _useCaloHit(conf().useCaloHit()),
    _debug(conf().debug())
  {
    produces<KalSeedPtrCollection>();
    if(_makeHists) {
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir("resolving");
      _hists.fitcon = tfdir.make<TH1F>("fitcon", "Track p(#chi^2)", 200,  0.,  1.);
    }
  }

  //--------------------------------------------------------------------------------------------------------
  // Decide which track is better
  int TrackResolution::resolve(const KalSeed* k_1, const KalSeed* k_2) {
    if(!k_1 && !k_2) return kNeither;
    if(!k_1) return kSecond;
    if(!k_2) return kFirst;

    const float fitcon_1(k_1->fitConsistency()), fitcon_2(k_2->fitConsistency());
    if(_debug > 1) printf("    Comparing two tracks: Calo cluster: k1 = %o, k2 = %o; Fit quality: k1 = %.3g, k2 = %.3g\n",
                          k_1->hasCaloCluster(), k_2->hasCaloCluster(), fitcon_1, fitcon_2);

    // If requested, prioritize tracks that include a calorimeter cluster
    if(_useCaloHit ) {
      if( k_1->hasCaloCluster() && !k_2->hasCaloCluster()) return kFirst;
      if(!k_1->hasCaloCluster() &&  k_2->hasCaloCluster()) return kSecond;
    }

    // Compare the fit quality
    if(fitcon_1 > fitcon_2) return kFirst;
    return kSecond;
  }

  //--------------------------------------------------------------------------------------------------------
  void TrackResolution::produce(art::Event& event ) {
    // create output
    std::unique_ptr<KalSeedPtrCollection> trkcol(new KalSeedPtrCollection());

    // get the KalSeedPtrs
    auto handle = event.getValidHandle<KalSeedPtrCollections>(_seedTag);
    const auto& clusters = *handle;
    if(_debug > 0) printf("[TrackResolution::%s::%s] Collection of %zu track clusters retrieved\n", __func__, moduleDescription().moduleLabel().c_str(), clusters.size());

    // Go through the tracks and check for alternate track fit hypotheses using the same helix seed
    for (const auto& cluster : clusters) {
      if(cluster.empty()) continue;
      if(_debug > 1) printf("  Checking cluster of %zu tracks\n", cluster.size());
      const KalSeedPtr* seedPtr = &(cluster.front());
      if(_makeHists) {
        _hists.fitcon->Fill((*seedPtr)->fitConsistency());
      }

      // Resolve the overlap for each track in the cluster
      const size_t ntrks = cluster.size();
      for(size_t itrk = 1; itrk < ntrks; ++itrk) {
        const KalSeedPtr testPtr = cluster.at(itrk);
        const KalSeed* k_1 = &(**seedPtr);
        const KalSeed* k_2 = &(*testPtr);
        const int result = resolve(k_1, k_2);
        if(_debug > 1) printf("  Resolving against track %zu has status %i\n", itrk, result);
        if(result == kSecond) seedPtr = &testPtr;
        else if(result == kNeither) seedPtr = nullptr;
        if(_makeHists) {
          _hists.fitcon->Fill((*seedPtr)->fitConsistency());
        }
      }

      if(seedPtr) {
        if(_debug > 1) {
          // Approximate the fit trajectory by just taking the first tracker segment pZ sign
          int traj = 0;
          float mom = 0.f;
          for(auto inter : (*seedPtr)->intersections()) {
            if(inter.surfaceId() == SurfaceIdDetail::TT_Mid ||
               inter.surfaceId() == SurfaceIdDetail::TT_Front ||
               inter.surfaceId() == SurfaceIdDetail::TT_Back) {
              traj = (inter.momentum3().z() < 0.) ? -1 : 1;
              mom = inter.mom();
              break;
            }
          }
          printf("  Final track resolution: PDG = %5i, Trajectory = %2i, p = %5.1f\n", (*seedPtr)->particle(), traj, mom);
        }
        trkcol->push_back(*seedPtr);
      } else {
        if(_debug > 0) printf("[TrackResolution::%s::%s] A track cluster did not resolve to any track\n", __func__, moduleDescription().moduleLabel().c_str());
      }
    }

    if(_debug > 0) printf("[TrackResolution::%s::%s] Outputting a track collection of size %zu\n", __func__, moduleDescription().moduleLabel().c_str(), trkcol->size());

    // put the output products into the event
    event.put(std::move(trkcol));
  }
}// mu2e

DEFINE_ART_MODULE(mu2e::TrackResolution)
