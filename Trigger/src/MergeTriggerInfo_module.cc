//
//  Merge the TriggerInfo objects generated within the trigger-path
//  author: G. Pezzullo
//
#include "fhiclcpp/types/Atom.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
// mu2e data products
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloTrigSeed.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
// C++
#include <vector>
#include <memory>
#include <iostream>
#include <forward_list>
#include <string>


namespace mu2e {
  class MergeTriggerInfo : public art::EDProducer {
  public:
    using  Name    = fhicl::Name;
    using  Comment = fhicl::Comment;
    struct Config {
      fhicl::Atom<int> debug     { Name("debugLevel"), Comment("Debug Level"), 0};
      fhicl::Atom<int> doDeepCopy{ Name("doDeepCopy"), Comment("Produce cloned object collections"), 0};
    };

    using        Parameters = art::EDProducer::Table<Config>;
    explicit     MergeTriggerInfo(const Parameters& conf);
    void         produce(art::Event& evt) override;
  private:
    int          _debug;
    int          _doDeepCopy;
  };

  MergeTriggerInfo::MergeTriggerInfo(const Parameters& config) :
    art::EDProducer{config},
    _debug     (config().debug()),
    _doDeepCopy(config().doDeepCopy())
  {
    produces<TriggerInfoCollection>();
    if (_doDeepCopy == 1){
      produces<KalSeedCollection        >();
      produces<HelixSeedCollection      >();
      produces<TimeClusterCollection    >();
      produces<CaloClusterCollection    >();
      produces<CaloTrigSeedCollection   >();
      produces<CosmicTrackSeedCollection>();
    }
  }

  void MergeTriggerInfo::produce(art::Event& event) {
    std::unique_ptr<TriggerInfoCollection>     tiCol(new TriggerInfoCollection);
    std::unique_ptr<KalSeedCollection>         ksCol((_doDeepCopy) ? new KalSeedCollection         : nullptr);
    std::unique_ptr<HelixSeedCollection>       hsCol((_doDeepCopy) ? new HelixSeedCollection       : nullptr);
    std::unique_ptr<TimeClusterCollection>     tcCol((_doDeepCopy) ? new TimeClusterCollection     : nullptr);
    std::unique_ptr<CaloClusterCollection>     ccCol((_doDeepCopy) ? new CaloClusterCollection     : nullptr);
    std::unique_ptr<CaloTrigSeedCollection>    ctCol((_doDeepCopy) ? new CaloTrigSeedCollection    : nullptr);
    std::unique_ptr<CosmicTrackSeedCollection> csCol((_doDeepCopy) ? new CosmicTrackSeedCollection : nullptr);

    std::vector<art::Handle<TriggerInfo> > list_of_triggerInfo = event.getMany<TriggerInfo>();

    if(_debug > 0){
      std::cout << "["<<moduleDescription().moduleLabel() << "] number of TriggerInfo found in the is: "<< list_of_triggerInfo.size() << std::endl;
    }

    for (auto & trigInfoH: list_of_triggerInfo){
      TriggerInfo trigInfo(*trigInfoH.product());
      if(_debug > 0){
        std::cout << "["<<moduleDescription().moduleLabel() << "] helices, tracks: "<< trigInfo.tracks().size() << " " << trigInfo.helixes().size() << std::endl;
      }
      if (_doDeepCopy == 1){
        for(auto ptr : trigInfo.tracks       ()) ksCol->push_back(*(ptr.get()));
        for(auto ptr : trigInfo.helixes      ()) hsCol->push_back(*(ptr.get()));
        for(auto ptr : trigInfo.hitClusters  ()) tcCol->push_back(*(ptr.get()));
        for(auto ptr : trigInfo.caloClusters ()) ccCol->push_back(*(ptr.get()));
        for(auto ptr : trigInfo.caloTrigSeeds()) ctCol->push_back(*(ptr.get()));
        for(auto ptr : trigInfo.cosmics      ()) csCol->push_back(*(ptr.get()));
      }
      tiCol->push_back(trigInfo);
    }

    event.put(std::move(tiCol));
    if (_doDeepCopy == 1){
      event.put(std::move(ksCol));
      event.put(std::move(hsCol));
      event.put(std::move(tcCol));
      event.put(std::move(ccCol));
      event.put(std::move(ctCol));
      event.put(std::move(csCol));
    }
  }
}
DEFINE_ART_MODULE(mu2e::MergeTriggerInfo)
