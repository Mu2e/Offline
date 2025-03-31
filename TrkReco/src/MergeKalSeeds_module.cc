//
//  Merge KalSeeds from different reconstruction paths.
//  This module produces an indirect (Ptr) collection, not a deep copy, to keep from break provenance or MC truth associations
//  Note that, when persisting merged collections created by this module, the KalSeed object collections merged here MUST ALSO BE PERSISTED in
//  order for the references to not be broken.
//  original author: D. Brown (LBNL) 2024
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalDelegatedParameter.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "Offline/TrkReco/inc/KalSeedSelector.hh"
// mu2e data products
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
// C++
#include <vector>

namespace mu2e {
  // sorting function
  class SortKalSeeds    {
    public:
      SortKalSeeds(KalSeedSelector const& sel): sel_(sel) {}
      bool operator()(KalSeedPtr const& a, KalSeedPtr const& b) const { return sel_.isBetter(*b,*a); }
    private:
      KalSeedSelector const& sel_;
  };


  class MergeKalSeeds : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<int> debug{ Name("debugLevel"), Comment("Debug Level"), 0};
        fhicl::Sequence<art::InputTag> seedCollTags {Name("KalSeedCollections"), Comment("KalSeed collections tomerge") };
        fhicl::Atom<bool> selbest{ Name("SelectBest"), Comment("Select best candidate"), false};
        fhicl::OptionalDelegatedParameter selector{Name("Selector"), Comment("Selector parameters")};
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit MergeKalSeeds(const Parameters& config);
      void produce(art::Event& evt) override;
    private:
      int debug_;
      bool selbest_;
      std::vector<art::InputTag> seedcolltags_;
      std::unique_ptr<KalSeedSelector> selector_;
  };

  MergeKalSeeds::MergeKalSeeds(const Parameters& config) : art::EDProducer{config},
    debug_(config().debug()),
    selbest_(config().selbest()) {
      for(auto const& seedcolltag :config().seedCollTags()){
        mayConsume<KalSeedCollection>    (seedcolltag);
        seedcolltags_.push_back(seedcolltag);
        if(debug_ > 0) std::cout << "Merging KalSeeds from " << seedcolltag << std::endl;
      }
      // if requested, instantiate the selector
      const auto pset = config().selector.get_if_present<fhicl::ParameterSet>();
      if(pset.has_value()) selector_ = art::make_tool<KalSeedSelector>(*pset);
      produces<KalSeedPtrCollection> ();
    }

  void MergeKalSeeds::produce(art::Event& event) {
    // create output
    std::unique_ptr<KalSeedPtrCollection> mkseeds(new KalSeedPtrCollection);
    // loop over input KalSeed collections
    for(auto const& seedcolltag : seedcolltags_) {
      auto const& ksch = event.getHandle<KalSeedCollection>(seedcolltag);
      if(ksch.isValid()){
        auto const& ksc = *ksch;
        if(debug_ > 2) std::cout << seedcolltag << " Has " << ksc.size() << " KalSeeds" << std::endl;
        for(size_t ikseed = 0; ikseed <ksc.size(); ++ikseed){
          auto const& kseed = ksc[ikseed];
          if( !selector_ ) {
            mkseeds->emplace_back(ksch,ikseed);
          } else if ( selector_->select(kseed)){
            if(!selbest_ || mkseeds->size() == 0)
              mkseeds->emplace_back(ksch,ikseed);
            else if(selbest_ && selector_->isBetter(*mkseeds->back(),kseed)){
              mkseeds->back() = art::Ptr<KalSeed>(ksch,ikseed);
            }
          }
        }
      } else if (debug_ > 2){
        std::cout << "No collection found for " << seedcolltag << std::endl;
      }
    }
    // if we aren't selecting the best but a selector was provided, sort the ptrs
    if(selector_ && !selbest_){
      SortKalSeeds sortks(*selector_);
      std::sort(mkseeds->begin(),mkseeds->end(),sortks);
    }
    if(debug_ > 1) std::cout << "Merged " << mkseeds->size() << " KalSeeds" << std::endl;
    event.put(std::move(mkseeds));
  }
}
DEFINE_ART_MODULE(mu2e::MergeKalSeeds)
