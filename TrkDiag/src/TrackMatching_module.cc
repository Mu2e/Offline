//
// Cluster alternate fits of tracks and tracks with significant overlap
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
      fhicl::Atom<double> timeWindow{Name("timeWindow"), Comment("Hit time window in matching, if not using indices"), 50.};
      fhicl::Atom<bool> makeHists{Name("makeHistograms"), Comment("Make debug histograms"), false};
      fhicl::Atom<int> debug{Name("debugLevel"), Comment("Debug printout level"), 0};
    };

    using Parameters = art::EDProducer::Table<Config>;
    TrackMatching(const Parameters& conf);

  private:
    void produce(art::Event& event) override;
    bool match(const KalSeed* k_1, const KalSeed* k_2);
    void add_tracks(std::vector<std::vector<const KalSeed*>>& clusters,
                    std::map<const KalSeed*, size_t>& cluster_index,
                    const KalSeed* itrkptr, const KalSeed* jtrkptr);
    void create_clusters(std::vector<const KalSeedCollection*>& trackColls,
                         std::vector<art::ValidHandle<KalSeedCollection>>& handles,
                         std::map<const KalSeed*, KalSeedPtr>& ptr_map,
                         std::vector<std::vector<const KalSeed*>>& clusters,
                         std::map<const KalSeed*, size_t>& cluster_index);

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
    _timeWindow(conf().timeWindow()),
    _makeHists(conf().makeHists()),
    _debug(conf().debug())
  {
    produces<KalSeedPtrCollections>();
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
  bool TrackMatching::match(const KalSeed* k_1, const KalSeed* k_2) {
    if(!k_1 || !k_2)
      throw cet::exception("RECO") << "mu2e::TrackMatching::" << __func__ << ": Null input track seeds!";

    // Retrieve the hit lists
    const auto hits_1 = k_1->hits();
    const auto hits_2 = k_2->hits();
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
    const bool pass = fraction > _threshold;
    if(_debug > 2) printf("  Track overlap: %2i hits, %5.3f fraction --> pass = %o\n", overlap, fraction, pass);
    if(_makeHists) {
      _hists.overlap->Fill(fraction);
    }
    return pass;
  }

  //--------------------------------------------------------------------------------------------------------
  // Add a track to the clusters, handling merging if needed
  void TrackMatching::add_tracks(std::vector<std::vector<const KalSeed*>>& clusters,
                                std::map<const KalSeed*, size_t>& cluster_index,
                                const KalSeed* itrkptr, const KalSeed* jtrkptr) {
    const bool i_clustered = cluster_index.count(itrkptr);
    const bool j_clustered = cluster_index.count(jtrkptr);
    // Check if the tracks have already been clustered
    if(i_clustered && j_clustered) { // both are clustered --> merge their clusters
      const size_t i_index = cluster_index[itrkptr];
      const size_t j_index = cluster_index[jtrkptr];
      auto& i_cluster = clusters[i_index];
      auto& j_cluster = clusters[j_index];
      if(_debug > 2) printf("  Merging track clusters: %zu with %zu tracks and %zu with %zu tracks\n",
                            i_index, i_cluster.size(), j_index, j_cluster.size());

      // Add all of the jtrk cluster to the itrk cluster, then erase the jtrk entry
      for(auto ptr : j_cluster) {
        cluster_index[ptr] = i_index;
        i_cluster.push_back(ptr);
      }
      clusters[j_index] = {}; // replace with empty list to preserve the indexing
    } else if(i_clustered) { // itrk has been clustered
      const size_t index = cluster_index[itrkptr];
      clusters[index].push_back(jtrkptr);
      cluster_index[jtrkptr] = index;
    } else if(cluster_index.count(jtrkptr)) { // jtrk has been clustered
      const size_t index = cluster_index[jtrkptr];
      clusters[index].push_back(itrkptr);
      cluster_index[itrkptr] = index;
    } else { // Neither track has been clustered
      const size_t index = clusters.size();
      clusters.push_back({itrkptr, jtrkptr});
      cluster_index[itrkptr] = index;
      cluster_index[jtrkptr] = index;
    }
  }

  //--------------------------------------------------------------------------------------------------------
  // Create track clusters from input track collections
  void TrackMatching::create_clusters(std::vector<const KalSeedCollection*>& trackColls,
                                      std::vector<art::ValidHandle<KalSeedCollection>>& handles,
                                      std::map<const KalSeed*, KalSeedPtr>& ptr_map,
                                      std::vector<std::vector<const KalSeed*>>& clusters,
                                      std::map<const KalSeed*, size_t>& cluster_index) {

    const size_t ncolls = trackColls.size();
    const size_t max_coll = (ncolls < 1) ? 0 : ncolls - 1; // no need to check the last collection as others checked against it
    if(_debug > 2) printf("[TrackMatching::%s::%s] Inspecting track inputs from %zu collections\n",
                          __func__, moduleDescription().moduleLabel().c_str(), ncolls);

    for(size_t i = 0; i < max_coll; ++i) {
      const KalSeedCollection* itrks = trackColls[i];
      const size_t n_itrks = itrks->size();
      if(_debug > 2) printf("  Checking track collection %s (size = %zu)\n", _seedTags[i].encode().c_str(), n_itrks);

      // Check each track for overlaps
      for(size_t i_index = 0; i_index < n_itrks; ++i_index) {
        if(_debug > 3) printf("    Checking track %zu of collection %s\n", i_index, _seedTags[i].encode().c_str());
        const auto& itrk = itrks->at(i_index);
        const auto* itrkptr = &itrk; // pointer to the track
        if(_debug > 4) printf("    Retrieved the track\n");
        if(!ptr_map.count(itrkptr)) { // add the track to the map
          if(_debug > 3) printf("  --> Adding track %zu of %s to the Ptr map\n", i_index, _seedTags[i].encode().c_str());
          ptr_map[itrkptr] = KalSeedPtr(handles[i], i_index);
        }

        // Check against each track collection (including its own collection), ignoring earlier collections already checked
        for(size_t j = i; j < ncolls; ++j) {
          const KalSeedCollection* jtrks = trackColls[j];
          const size_t n_jtrks = jtrks->size();
          if(_debug > 2) printf("    Checking against track collection %s (size = %zu)\n", _seedTags[j].encode().c_str(), n_jtrks);

          // Check against each track
          for(size_t j_index = 0; j_index < n_jtrks; ++j_index) {
            if(_debug > 3) printf("      Checking track %zu of collection %s\n", j_index, _seedTags[j].encode().c_str());
            const auto& jtrk = jtrks->at(j_index);
            const auto* jtrkptr = &jtrk; // pointer to the track

            // Check if both tracks have already been clustered together by another track
            const bool i_clustered = cluster_index.count(itrkptr);
            const bool j_clustered = cluster_index.count(jtrkptr);
            if(i_clustered && j_clustered && cluster_index[itrkptr] == cluster_index[jtrkptr]) continue;

            // Ensure this is not the same track
            if(itrkptr == jtrkptr) {
              if(i != j && _debug > 0) printf("[TrackMatching::%s::%s] A track is in two input track collections! Collection %s:%zu and %s:%zu\n",
                                              __func__, moduleDescription().moduleLabel().c_str(),
                                              _seedTags[i].encode().c_str(), i_index,
                                              _seedTags[j].encode().c_str(), j_index);
              continue;
            }

            // If the tracks overlap, add the tracks to the match lists
            if(match(itrkptr, jtrkptr)) {
              if(_debug > 1) printf("[TrackMatching::%s::%s] Found a match!\n",
                                    __func__, moduleDescription().moduleLabel().c_str());
              if(!ptr_map.count(jtrkptr)) { // add the track to the map
                if(_debug > 3) printf("  --> Adding track %zu of %s to the Ptr map\n", j_index, _seedTags[j].encode().c_str());
                ptr_map[jtrkptr] = KalSeedPtr(handles[j], j_index);
              }

              // add the tracks to the clusters
              add_tracks(clusters, cluster_index, itrkptr, jtrkptr);
            }
          }
        }
      }
    }
  }

  //--------------------------------------------------------------------------------------------------------
  void TrackMatching::produce(art::Event& event ) {
    // Create the output
    std::unique_ptr<KalSeedPtrCollections> ptr_clusters(new KalSeedPtrCollections);

    // Get the track collections
    std::vector<const KalSeedCollection*> trackColls;
    std::vector<art::ValidHandle<KalSeedCollection>> handles;
    for(auto tag : _seedTags) {
      auto handle = event.getValidHandle<KalSeedCollection>(tag);
      handles.push_back(handle);
      trackColls.push_back(handle.product());
      if(_debug > 0) printf("[TrackMatching::%s::%s] Track collection %s has %zu entries\n",
                            __func__, moduleDescription().moduleLabel().c_str(),
                            tag.encode().c_str(), handle->size());
    }

    //
    // For each track in each collection, check for overlaps with each other track
    //

    std::map<const KalSeed*, KalSeedPtr> ptr_map; // map of track -> art::Ptr
    std::vector<std::vector<const KalSeed*>> clusters; // list of clusters of matched tracks
    std::map<const KalSeed*, size_t> cluster_index; // map of track -> index in the cluster list

    create_clusters(trackColls, handles, ptr_map, clusters, cluster_index);

    //
    // From the clustered track lists, create the output track cluster collections
    //

    // Create the matched clusters
    int nclusters = 0; // count the non-empty clusters
    for(auto& cluster : clusters) {
      if(cluster.empty()) continue; // skip emptied clusters due to cluster merging
      if(_debug > 0) printf("  Adding a track cluster with %zu tracks\n", cluster.size());
      if(_makeHists) _hists.cluster_sizes->Fill(cluster.size());
      ++nclusters;
      KalSeedPtrCollection ptrs;
      for(auto trk : cluster) ptrs.push_back(ptr_map[trk]);
      ptr_clusters->push_back(ptrs);
    }


    // Add non-overlapping tracks as single track clusters so all tracks are in the output
    int nsingle_clusters = 0; // count the isolated tracks
    for(auto entry : ptr_map) {
      if(cluster_index.count(entry.first)) continue;
      if(_debug > 0) printf("  Adding a single track cluster\n");
      ++nsingle_clusters;
      KalSeedPtrCollection ptrs;
      ptrs.push_back(entry.second);
      if(_makeHists) _hists.cluster_sizes->Fill(1);
      ptr_clusters->push_back(ptrs);
    }

    if(_debug > 0) printf("[TrackMatching::%s::%s] Found %i track clusters and %i single tracks\n",
                          __func__, moduleDescription().moduleLabel().c_str(), nclusters, nsingle_clusters);
    if(_makeHists) _hists.nclusters->Fill(nclusters + nsingle_clusters);

    // Put the output products into the event
    event.put(std::move(ptr_clusters));
  }
}// mu2e

DEFINE_ART_MODULE(mu2e::TrackMatching)
