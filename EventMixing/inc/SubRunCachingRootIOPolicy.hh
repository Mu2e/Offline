#ifndef EventMixing_inc_SubRunCachingRootIOPolicy_hh
#define EventMixing_inc_SubRunCachingRootIOPolicy_hh
//
// An art mixing I/O policy that reads SubRun- and Run-level mixed products once per secondary SubRun (Run)
// instead of once per primary event.
//
// art's MixHelper asks the I/O policy for the SubRun products of the secondary events of every primary event,
// and art::RootIOPolicy reads them from the file every time. A mixer that draws many consecutive events from
// one secondary SubRun, such as the cosmic ResamplingMixer, reads the same products again for every event; for
// the stage-1 PhysicalVolumeInfoMultiCollection that costs more than reading the event products.
//
// This policy forwards everything to art::RootIOPolicy. For a SubRun- or Run-level mix op it keeps the products
// of the last read and returns them while the requested entries are the same. Opening a secondary file clears
// the cache. The products are const and the same entries of the same file hold the same products, so the mix
// functions see exactly what a fresh read would give them.
//

#include <map>
#include <string>

#include "art/Framework/IO/ProductMix/MixIOPolicy.h"
#include "art/Framework/IO/ProductMix/MixOpBase.h"
#include "art_root_io/RootIOPolicy.h"

namespace mu2e {

  class SubRunCachingRootIOPolicy : public art::MixIOPolicy {
  public:
    art::EventAuxiliarySequence
    generateEventAuxiliarySequence(art::EntryNumberSequence const& seq) override {
      return root().generateEventAuxiliarySequence(seq);
    }
    bool fileOpen() const override { return root().fileOpen(); }
    std::size_t nEventsInFile() const override { return root().nEventsInFile(); }
    art::FileIndex const& fileIndex() const override { return root().fileIndex(); }
    cet::exempt_ptr<art::BranchIDLists const> branchIDLists() const override { return root().branchIDLists(); }

    void openAndReadMetaData(std::string fileName, art::MixOpList& mixOps) override {
      cache_.clear();
      root().openAndReadMetaData(std::move(fileName), mixOps);
    }

    art::SpecProdList readFromFile(art::MixOpBase const& mixOp, art::EntryNumberSequence const& seq) override {
      if(mixOp.branchType() == art::InEvent || seq.empty()) {
        return root().readFromFile(mixOp, seq);
      }
      auto& c = cache_[&mixOp];
      if(!c.valid || c.seq != seq) {
        c.products = root().readFromFile(mixOp, seq);
        c.seq = seq;
        c.valid = true;
      }
      return c.products;
    }

  private:
    struct Cached {
      bool valid{false};
      art::EntryNumberSequence seq;
      art::SpecProdList products;
    };

    art::RootIOPolicy root_;
    std::map<art::MixOpBase const*, Cached> cache_;

    // art::RootIOPolicy declares its overrides private: call them through the public base interface.
    art::MixIOPolicy& root() { return root_; }
    art::MixIOPolicy const& root() const { return root_; }
  };

}

#endif /* EventMixing_inc_SubRunCachingRootIOPolicy_hh */
