//
//  Create Ptr collections from a trigger info collection or single trigger info object
//  author: Michael MacKenzie, 2025
//

// Framework
#include "fhiclcpp/types/Atom.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// Offline
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"

// C++
#include <vector>
#include <string>

namespace mu2e {
  class TriggerInfoToCollections : public art::EDProducer {
  public:
   struct Config {
      using  Name    = fhicl::Name;
      using  Comment = fhicl::Comment;
      fhicl::Atom<std::string> tag   { Name("triggerInfoName") , Comment("Trigger info (collection) name")};
      fhicl::Atom<int>         debug { Name("debugLevel")      , Comment("Debug Level"), 0};
    };

    using        Parameters = art::EDProducer::Table<Config>;
    explicit     TriggerInfoToCollections(const Parameters& conf);
    void         produce(art::Event& evt) override;

    void         fillContainers(const TriggerInfo& info,
                                KalSeedPtrCollection&     trkCol,
                                HelixSeedPtrCollection&   hlxCol,
                                TimeClusterPtrCollection& tclCol);
  private:
    std::string  _tag;
    int          _debug;
  };

  TriggerInfoToCollections::TriggerInfoToCollections(const Parameters& config) :
    art::EDProducer{config}
    , _tag(config().tag())
    , _debug(config().debug())
  {
    produces<KalSeedPtrCollection    >();
    produces<HelixSeedPtrCollection  >();
    produces<TimeClusterPtrCollection>();
  }

  void TriggerInfoToCollections::fillContainers(const TriggerInfo& info,
                                                KalSeedPtrCollection&     trkCol,
                                                HelixSeedPtrCollection&   hlxCol,
                                                TimeClusterPtrCollection& tclCol) {
    for(auto ptr : info.tracks()) {
      trkCol.push_back(ptr);
    }
    for(auto ptr : info.helixes()) {
      hlxCol.push_back(ptr);
    }
    for(auto ptr : info.hitClusters()) {
      tclCol.push_back(ptr);
    }
  }

  void TriggerInfoToCollections::produce(art::Event& event) {
    std::unique_ptr<KalSeedPtrCollection    > trkCol(new KalSeedPtrCollection    );
    std::unique_ptr<HelixSeedPtrCollection  > hlxCol(new HelixSeedPtrCollection  );
    std::unique_ptr<TimeClusterPtrCollection> tclCol(new TimeClusterPtrCollection);
    art::Handle<TriggerInfoCollection> info_coll_handle; // output from MergeTriggerInfo
    art::Handle<TriggerInfo> info_handle; // output from a single trigger filter instance

    if(event.getByLabel(_tag, info_coll_handle)) { // assuming we're looking for a collection of TriggerInfo objects
      auto infos = info_coll_handle.product();

      // fill the output Ptr collections with the trigger info collections
      for(const auto& info : *infos) {
        fillContainers(info, *trkCol, *hlxCol, *tclCol);
      }

      // report information if requested
      if(_debug > 1 || (_debug > 0 && !trkCol->empty())) {
        printf("[TriggerInfoToCollections::%s::%s] From %2zu objects found %2zu tracks, %2zu helices, and %2zu time clusters\n",
               __func__, moduleDescription().moduleLabel().c_str(), infos->size(), trkCol->size(), hlxCol->size(), tclCol->size());
      }
    } else if(event.getByLabel(_tag, info_handle)) { // check for a single TriggerInfo object by this label
      fillContainers(*info_handle.product(), *trkCol, *hlxCol, *tclCol);

      // report information if requested
      if(_debug > 1 || (_debug > 0 && !trkCol->empty())) {
        printf("[TriggerInfoToCollections::%s::%s] From single object found %2zu tracks, %2zu helices, and %2zu time clusters\n",
               __func__, moduleDescription().moduleLabel().c_str(), trkCol->size(), hlxCol->size(), tclCol->size());
      }
    } else {
      if(_debug > 1) {
        printf("[TriggerInfoToCollections::%s::%s] No TriggerInfo collection found with tag %s\n",
               __func__, moduleDescription().moduleLabel().c_str(), _tag.c_str());
      }
    }

    // add the collections to the event
    event.put(std::move(trkCol));
    event.put(std::move(hlxCol));
    event.put(std::move(tclCol));
  }
}
DEFINE_ART_MODULE(mu2e::TriggerInfoToCollections)
