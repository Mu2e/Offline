//
//  Merge HelixSeeds from different pattern recognition paths.  Modeled
//  on MergeHelixFinder.  P. Murat, D. Brown (LBNL)
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
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
// utilities
#include "TrkReco/inc/TrkUtilities.hh"
// C++
#include <vector>
#include <memory>
#include <iostream>
#include <set>
#include <string>


namespace mu2e {
  class MergeHelices : public art::EDProducer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Atom<int> debug{ Name("debugLevel"),
	Comment("Debug Level"), 0};
      fhicl::Atom<bool> selectbest{ Name("SelectBest"),
	Comment("Select best overlapping helices for output"), true};
      fhicl::Atom<bool> usecalo{ Name("UseCalo"),
	Comment("Use CaloCluster info in comparison"), true};
      fhicl::Atom<unsigned> minnover{ Name("MinNHitOverlap"),
	Comment("Minimum number of common hits to consider helices to be 'the same'"), 10};
      fhicl::Atom<float> minoverfrac{ Name("MinHitOverlapFraction"),
	Comment("Minimum fraction of common hits to consider helices to be 'the same'"), 0.5};
      fhicl::Sequence<std::string> BadHitFlags { Name("BadHitFlags"),
	Comment("HelixHit flag bits to exclude from counting"),std::vector<std::string>{"Outlier"}}; 
      fhicl::Sequence<std::string> HelixFinders { Name("HelixFinders"),
	Comment("HelixSeed producers to merge")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit MergeHelices(const Parameters& conf);
    void produce(art::Event& evt) override;
  private:
    enum HelixComp{unique=-1,first=0,second=1};
    int _debug;
    unsigned _minnover;
    float _minoverfrac;
    bool _selectbest, _usecalo;
    std::vector<std::string> _hfs;
    StrawHitFlag _badhit;
    // helper functions
    typedef std::vector<StrawHitIndex> SHIV;
    HelixComp compareHelices(art::Event const& evt,
	HelixSeed const& h1, HelixSeed const& h2);
    void countHits(art::Event const& evt,
	HelixSeed const& h1, HelixSeed const& h2,
	unsigned& nh1, unsigned& nh2, unsigned& nover);
    unsigned countOverlaps(SHIV const& s1, SHIV const& s2);
  };

  MergeHelices::MergeHelices(const Parameters& config) : 
    art::EDProducer{config},
    _debug(config().debug()),
    _minnover(config().minnover()),
    _minoverfrac(config().minoverfrac()),
    _selectbest(config().selectbest()),
    _usecalo(config().usecalo()),
    _hfs(config().HelixFinders()),
    _badhit(config().BadHitFlags())
  {
    consumesMany<HelixSeedCollection>    ();
    produces<HelixSeedCollection>    ();
    produces<TimeClusterCollection>    ();
  }

  void MergeHelices::produce(art::Event& event) {
  // create output
    std::unique_ptr<HelixSeedCollection> mhels(new HelixSeedCollection);
    std::unique_ptr<TimeClusterCollection> tcs(new TimeClusterCollection);
    // needed for creating Ptrs
    auto TimeClusterCollectionPID = event.getProductID<TimeClusterCollection>();
    auto TimeClusterCollectionGetter = event.productGetter(TimeClusterCollectionPID);
    // loop over helix products and flatten the helix collections into a single collection
    std::set<const HelixSeed*> hseeds;
    for (auto const& hf : _hfs) {
      art::InputTag hsct(hf);
      auto hsch = event.getValidHandle<HelixSeedCollection>(hsct);
      auto const& hsc = *hsch;
      for(auto const&  hs : hsc) {
	hseeds.insert(&hs);
      }
    }
// now loop over all combinations
    for(auto ihel = hseeds.begin(); ihel != hseeds.end();) {
      auto jhel = ihel; jhel++;
      while( jhel != hseeds.end()){
	// compare the helice 
	auto hcomp = compareHelices(event, **ihel, **jhel);
	if(hcomp == unique) {
	  // both helices are unique: simply advance the iterator to keep both
	  jhel++;
	} else if(hcomp == first) {
	  // the first helix is 'better'; remove the second
	  jhel = hseeds.erase(jhel);
	} else if(hcomp == second) {
	  // the second helix is better; remove the first and restart the loop
	  ihel = hseeds.erase(ihel);
	  break;
	}
      }
      // only advance the outer loop if we exhausted the inner one
      if(jhel == hseeds.end())ihel++;
    }
    // now hseeds contains only pointers to unique and best helices: use them to 
    // create the output, which must be deep-copy, including the time cluster
    for(auto ihel = hseeds.begin(); ihel != hseeds.end(); ihel++){
      // copy the time cluster first
      tcs->push_back(*(*ihel)->timeCluster());
      // create a Ptr to the new TimeCluster
      auto tcptr = art::Ptr<TimeCluster>(TimeClusterCollectionPID,tcs->size()-1,TimeClusterCollectionGetter);
      // copy the Helix Seed and point its TimeCluster to the copy
      HelixSeed hs(**ihel);
      hs._timeCluster = tcptr;
      mhels->push_back(hs);
    }
    event.put(std::move(mhels));
    event.put(std::move(tcs));
  }

  MergeHelices::HelixComp MergeHelices::compareHelices(art::Event const& evt,
    HelixSeed const& h1, HelixSeed const& h2) {
    HelixComp retval(unique);
  // count the StrawHit overlap between the helices
    unsigned nh1, nh2, nover;
    countHits(evt,h1,h2, nh1, nh2, nover);
    unsigned minh = std::min(nh1, nh2);
    if(nover >= _minnover && nover/float(minh) > _minoverfrac) {
    // overlapping helices: decide which is best
    // this should be a neural net FIXME!
    // Pick the one with a CaloCluster first
      if(h1.caloCluster().isNonnull() && h2.caloCluster().isNull())
	retval = first;
      else if( h2.caloCluster().isNonnull() && h1.caloCluster().isNull())
	retval = first;
	// then compare active StrawHit counts
      else if(nh1 > nh2)
	retval = first;
      else if(nh2 > nh1)
	retval = second;
	// finally compare chisquared: sum xy and fz for now
      else if(h1.helix().chi2dXY() + h1.helix().chi2dZPhi()  < 
	  h2.helix().chi2dXY() + h2.helix().chi2dZPhi() )
	retval = first;
      else
	retval = second;
    }
    return retval;
  }

  void MergeHelices::countHits(art::Event const& evt,
    HelixSeed const& h1, HelixSeed const& h2,
    unsigned& nh1, unsigned& nh2, unsigned& nover) {
    nh1 = nh2 = nover = 0;
    SHIV shiv1, shiv2;
    for(size_t ihit=0;ihit < h1.hits().size(); ihit++){
      auto const& hh = h1.hits()[ihit];
      if(!hh.flag().hasAnyProperty(_badhit))
	h1.hits().fillStrawHitIndices(evt,ihit,shiv1);
    }
    for(size_t ihit=0;ihit < h2.hits().size(); ihit++){
      auto const& hh = h2.hits()[ihit];
      if(!hh.flag().hasAnyProperty(_badhit))
	h2.hits().fillStrawHitIndices(evt,ihit,shiv2);
    }
    nh1 = shiv1.size();
    nh2 = shiv2.size();
    nover = countOverlaps(shiv1,shiv2);
  }
  
  unsigned MergeHelices::countOverlaps(SHIV const& shiv1, SHIV const& shiv2) {
    unsigned retval(0);
    for(auto shi1 : shiv1)
      for(auto shi2 : shiv2)
	if(shi1 == shi2)retval++;
    return retval;
  }
}
DEFINE_ART_MODULE(mu2e::MergeHelices)
