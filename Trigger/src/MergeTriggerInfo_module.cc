//
//  Merge the TriggerInfo objects generated within the trigger-path
//  author: G. Pezzullo
//
#include "fhiclcpp/types/Atom.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
// mu2e data products
#include "RecoDataProducts/inc/TriggerInfo.hh"
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
      fhicl::Atom<int> debug{ Name("debugLevel"),
	  Comment("Debug Level"), 0};
    };

    using        Parameters = art::EDProducer::Table<Config>;
    explicit     MergeTriggerInfo(const Parameters& conf);
    void         produce(art::Event& evt) override;
  private:
    int          _debug;
  };

  MergeTriggerInfo::MergeTriggerInfo(const Parameters& config) : 
    art::EDProducer{config},
    _debug   (config().debug())
  {
    produces<TriggerInfoCollection>    ();
  }

  void MergeTriggerInfo::produce(art::Event& event) {
    std::unique_ptr<TriggerInfoCollection> tiCol(new TriggerInfoCollection);

    // create the selector
    art::Selector selector(art::ProductInstanceNameSelector("") &&
			   art::ProcessNameSelector("*")); 
    //			   art::ModuleLabelSelector(""));// && 
    //                           art::ProcessNameSelector("*"));
    std::vector<art::Handle<TriggerInfo> > list_of_triggerInfo; 
    event.getMany<TriggerInfo>(selector);

    if(_debug > 0){
      std::cout << "["<<moduleDescription().moduleLabel() << "] number of TriggerInfo found in the is: "<< list_of_triggerInfo.size() << std::endl;
    }
    
    for (auto & trigInfoH: list_of_triggerInfo){
      TriggerInfo trigInfo(*trigInfoH.product());
      tiCol->push_back(trigInfo);
    }
    
    event.put(std::move(tiCol));
  }
}
DEFINE_ART_MODULE(mu2e::MergeTriggerInfo)
