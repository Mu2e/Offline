//
// Cluster alternate fits of tracks and tracks with significant overlap
//
// Original author M. MacKenzie (2025)
//

// framework
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art_root_io/TFileService.h"

// data
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/TrkStrawHitSeed.hh"

// ROOT
#include "TFile.h"
#include "TH1.h"

// C++
#include <iostream>
#include <map>
#include <vector>

namespace mu2e
{

  class TrackMatching : public art::EDProducer
  {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Sequence<art::InputTag> seedTags{Name("TrackCollections"), Comment("Input tags for track collections")};
      fhicl::Atom<double> threshold{Name("overlapThreshold"), Comment("Fractional threshold in hits to match two tracks"), 0.5};
      fhicl::Atom<bool> useIndices{Name("useIndices"), Comment("Use digi indices when matching, otherwise use straw ID + time"), true};
      fhicl::OptionalAtom<double> timeWindow{Name("timeWindow"), Comment("Hit time window in matching, if not using indices -- default = 50 ns")};
      fhicl::Atom<bool> makeHists{Name("makeHistograms"), Comment("Make debug histograms"), false};
      fhicl::Atom<int> debug{Name("debugLevel"), Comment("Debug printout level"), 0};
    };

    using Parameters = art::EDProducer::Table<Config>;
    TrackMatching(const Parameters& conf);

  private:
    void produce(art::Event& event) override;
    bool match(const KalSeedPtr& k_1, const KalSeedPtr& k_2);
    void merge_clusters(KalSeedCluster& cluster, KalSeedCluster& cluster_j);
    void create_clusters(std::vector<KalSeedCluster>& clusters);

    std::vector<art::InputTag> _seedTags;
    double _threshold;
    bool _useIndices;
    double _timeWindow;
    bool _makeHists;
    int _debug;

    struct Hist_t {
      TH1* overlap;
      TH1* nclusters;
      TH1* cluster_sizes;
    };
    Hist_t _hists;
  };

  //--------------------------------------------------------------------------------------------------------
  TrackMatching::TrackMatching(const Parameters& conf) :
    art::EDProducer{conf},
    _seedTags(conf().seedTags()),
    _threshold(conf().threshold()),
    _useIndices(conf().useIndices()),
    _makeHists(conf().makeHists()),
    _debug(conf().debug())
  {
    produces<KalSeedClusterCollection>();
    if(!conf().timeWindow(_timeWindow)) _timeWindow = 50.; //default to a 50 ns window

    if(_makeHists) {
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir("matching");
      _hists.overlap = tfdir.make<TH1F>("overlap", "Track hit overlap", 110,  0.,  1.1);
      _hists.nclusters = tfdir.make<TH1F>("nclusters", "N(track clusters) / event", 20, 0, 20);
      _hists.cluster_sizes = tfdir.make<TH1F>("cluster_sizes", "N(tracks) / cluster", 20, 0, 20);
    }
  }

  //--------------------------------------------------------------------------------------------------------
  // Perform the matching
  bool TrackMatching::match(const KalSeedPtr& k_1, const KalSeedPtr& k_2) {
    if(_debug > 3) printf("  Starting track overlap check\n");
    if(k_1.isNull() || k_2.isNull())
      throw cet::exception("RECO") << "mu2e::TrackMatching::" << __func__ << ": Null input track seeds!";

    // Retrieve the hit lists
    const auto& hits_1 = k_1->hits();
    const auto& hits_2 = k_2->hits();
    if(hits_1.empty() || hits_2.empty()) return false;

    // Count the number of overlapping hits
    int overlap(0);
    for(auto hit_1 : hits_1) {
      for(auto hit_2 : hits_2) {
        if(_useIndices) { // match using the digi index
          if(hit_1._index == hit_2._index) ++overlap;
        } else {
          if(hit_1._sid == hit_2._sid && std::fabs(hit_1.hitTime() - hit_2.hitTime()) < _timeWindow) ++overlap;
        }
      }
    }
    const double fraction = overlap * 2. / (hits_1.size() + hits_2.size());
    const bool pass = fraction >= _threshold;
    if(_debug > 2) printf("  Track overlap: %2i hits, %5.3f fraction --> pass = %o\n", overlap, fraction, pass);
    if(_makeHists) {
      _hists.overlap->Fill(fraction);
    }
    return pass;
  }

  //--------------------------------------------------------------------------------------------------------
  // Merge two track clusters into the first cluster
  void TrackMatching::merge_clusters(KalSeedCluster& cluster, KalSeedCluster& cluster_j) {
    for(const auto& ptr : cluster_j) cluster.push_back(ptr); // merge into the first cluster
    cluster_j = {}; // empty the second cluster
  }

  //--------------------------------------------------------------------------------------------------------
  // Create track clusters from input track collections
  void TrackMatching::create_clusters(std::vector<KalSeedCluster>& clusters) {

    const size_t ntrks = clusters.size();
    if(_debug > 2) printf("[TrackMatching::%s::%s] Inspecting %zu tracks\n",
                          __func__, moduleDescription().moduleLabel().c_str(), ntrks);

    // Loop through the clusters, merging overlapping clusters
    for(size_t index = 0; index < ntrks; ++index) {
      KalSeedCluster& cluster = clusters[index];
      if(cluster.empty()) continue; // already merged into an earlier cluster
      if(cluster.size() > 1) // input clusters should always be 1 track, as matches are clustered into lower index clusters
        throw cet::exception("RECO") << "mu2e::TrackMatching::" << __func__ << ": Input track cluster has already been clustered! Size = " << cluster.size();
      if(_debug > 2) printf("  Checking track cluster %zu for overlaps\n", index);

      // Check each following cluster (of size 1) for overlapping tracks
      for(size_t jtrk = index + 1; jtrk < ntrks; ++jtrk) {
        KalSeedCluster& cluster_j = clusters[jtrk];
        if(cluster_j.empty()) continue;
        if(cluster_j.size() > 1) // input clusters should always be 1 track, as matches are clustered into lower index clusters
          throw cet::exception("RECO") << "mu2e::TrackMatching::" << __func__ << ": Track cluster matching against has already been clustered! Size = " << cluster_j.size();

        // Get the tracks to compare
        KalSeedPtr& k_i = cluster.front();
        KalSeedPtr& k_j = cluster_j.front();
        if(k_i.isNull() || k_j.isNull())
          throw cet::exception("RECO") << "mu2e::TrackMatching::" << __func__ << ": Bad track in clusters";

        if(_debug > 2) printf("    Checking against track cluster %zu for overlaps\n", jtrk);

        // If the tracks overlap, add the tracks to the match lists
        if(match(k_i, k_j)) {
          if(_debug > 1) printf("[TrackMatching::%s::%s] Found a match! Track %zu <--> %zu\n",
                                __func__, moduleDescription().moduleLabel().c_str(), index, jtrk);

          // merge the clusters
          merge_clusters(cluster, cluster_j);
        }
      }
    }
  }

  //--------------------------------------------------------------------------------------------------------
  void TrackMatching::produce(art::Event& event ) {

    // Create the output collection
    std::unique_ptr<KalSeedClusterCollection> out_clusters(new KalSeedClusterCollection);

    // Get the track collections
    std::vector<art::ValidHandle<KalSeedCollection>> handles;
    size_t ntrks = 0; // track how many tracks should be clustered
    for(auto tag : _seedTags) {
      auto handle = event.getValidHandle<KalSeedCollection>(tag);
      handles.push_back(handle);
      ntrks += handle->size();
      if(_debug > 0) printf("[TrackMatching::%s::%s] Track collection %s has %zu entries\n",
                            __func__, moduleDescription().moduleLabel().c_str(),
                            tag.encode().c_str(), handle->size());
    }

    // Start by creating a track cluster for each track, then merge collections as overlaps are found

    std::vector<KalSeedCluster> clusters;
    for(auto& handle : handles) {
      const size_t ntrks = handle->size();
      for(size_t itrk = 0; itrk < ntrks; ++itrk) {
        KalSeedPtr ptr(handle, itrk);
        clusters.push_back(KalSeedCluster({ptr}));
      }
    }

    // Perform the track matching and clustering

    create_clusters(clusters);


    // Add the clusters to the output

    size_t ntrks_out = 0;
    int nclusters(0), nsingle_clusters(0);
    for(auto& cluster : clusters) {
      if(cluster.empty()) continue; // skip emptied clusters from merging
      out_clusters->push_back(cluster);
      size_t ntrks_cl = cluster.size();
      ntrks_out += ntrks_cl;
      if(ntrks_cl == 1) ++nsingle_clusters;
      else              ++nclusters;
    }

    // Check that the right number of tracks were clustered
    if(ntrks != ntrks_out)
      throw cet::exception("RECO") << "mu2e::TrackMatching::" << __func__ << ": Number of clustered tracks doesn't match input! In = " << ntrks << " out = " << ntrks_out;

    if(_debug > 0) printf("[TrackMatching::%s::%s] Found %i track clusters and %i single tracks from %zu input tracks\n",
                          __func__, moduleDescription().moduleLabel().c_str(), nclusters, nsingle_clusters, ntrks);
    if(_makeHists) _hists.nclusters->Fill(nclusters + nsingle_clusters);

    // Put the output products into the event
    event.put(std::move(out_clusters));
  }
}// mu2e

DEFINE_ART_MODULE(mu2e::TrackMatching)
