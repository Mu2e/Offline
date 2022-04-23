//
//  Merge the IntensityInfo objects generated within the trigger-path
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
#include "Offline/RecoDataProducts/inc/IntensityInfo.hh"
// C++
#include <vector>
#include <memory>
#include <iostream>
#include <forward_list>
#include <string>


namespace mu2e {
  class MergeIntensityInfo : public art::EDProducer {
  public:
    using  Name    = fhicl::Name;
    using  Comment = fhicl::Comment;
    struct Config {
      fhicl::Atom<int> debug{ Name("debugLevel"),
	  Comment("Debug Level"), 0};
    };

    using        Parameters = art::EDProducer::Table<Config>;
    explicit     MergeIntensityInfo(const Parameters& conf);
    void         produce(art::Event& evt) override;
  private:
    int          _debug;
  };

  MergeIntensityInfo::MergeIntensityInfo(const Parameters& config) : 
    art::EDProducer{config},
    _debug   (config().debug())
  {
    produces<IntensityInfo>();
  }

  void MergeIntensityInfo::produce(art::Event& event) {
    std::unique_ptr<IntensityInfo> intInfo(new IntensityInfo);

    // create the selector
    art::Selector selector(art::ProductInstanceNameSelector("") &&
			   art::ProcessNameSelector("*")); 
    std::vector<art::Handle<IntensityInfo> > list_of_intensityInfo; 
    event.getMany<IntensityInfo>(selector);

    if(_debug > 0){
      std::cout << "["<<moduleDescription().moduleLabel() << "] number of IntensityInfo found in the is: "<< list_of_intensityInfo.size() << std::endl;
    }
    
    for (auto & evtInfoH: list_of_intensityInfo){
      IntensityInfo evtInfo(*evtInfoH.product());
      
      //loop over the data-members and store the values different from 0
      if (evtInfo.nTrackerHits() != 0){
	if (intInfo->nTrackerHits() != 0) std::cout << "["<<moduleDescription().moduleLabel() << "] nTrackerHits was already != 0 : "<< intInfo->nTrackerHits() << std::endl;
	intInfo->setNTrackerHits(evtInfo.nTrackerHits());
      }
      if (evtInfo.nCaloHits() != 0){
	if (intInfo->nCaloHits() != 0) std::cout << "["<<moduleDescription().moduleLabel() << "] nCaloHits was already != 0 : "<< intInfo->nCaloHits() << std::endl;
	intInfo->setNCaloHits(evtInfo.nCaloHits());
      }
      if (evtInfo.nProtonTCs() != 0){
	if (intInfo->nProtonTCs() != 0) std::cout << "["<<moduleDescription().moduleLabel() << "] nProtonTCs was already != 0 : "<< intInfo->nProtonTCs() << std::endl;
	intInfo->setNProtonTCs(evtInfo.nProtonTCs());
      }
      if (evtInfo.caloEnergy() != 0){
	if (intInfo->caloEnergy() != 0) std::cout << "["<<moduleDescription().moduleLabel() << "] caloEnergy was already != 0 : "<< intInfo->caloEnergy() << std::endl;
	intInfo->setCaloEnergy(evtInfo.caloEnergy());
      }
      if (evtInfo.testVariable() != 0){
	if (intInfo->testVariable() != 0) std::cout << "["<<moduleDescription().moduleLabel() << "] testVariable was already != 0 : "<< intInfo->testVariable() << std::endl;
	intInfo->setTestVariable(evtInfo.testVariable());
      }
    }
    
    event.put(std::move(intInfo));
  }
}
DEFINE_ART_MODULE(mu2e::MergeIntensityInfo)
