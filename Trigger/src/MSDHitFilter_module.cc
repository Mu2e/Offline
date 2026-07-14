//
//  Filter for mobile sync detector hits
//  Michael MacKenzie, 2026
//

// Framework
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"

// Offline
#include "Offline/RecoDataProducts/inc/MSDHit.hh"

// C++
#include <string>

// TRACE
#include "TRACE/tracemf.h"
#define TRACE_NAME "MSDHitFilter"

namespace mu2e {
  class MSDHitFilter : public art::EDFilter {
  public:

    // Main configuration
    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<std::string>  tag    {Name("tag")    , Comment("Collection tag")};
      fhicl::Atom<int>          minHits{Name("minHits"), Comment("Minimum number of hits")};
    };

    using Parameters = art::EDFilter::Table<Config>;

    explicit MSDHitFilter(const Parameters& config);


  private:
    bool filter  (art::Event& event) override;
    bool endRun  (art::Run&   run  ) override;
    bool goodHit (const MSDHit& hit);

    // Inputs
    std::string              _tag;
    int                      _minHits = -1;

    // Data
    unsigned long            _nevt, _npass;

  };

  //-----------------------------------------------------------------------------
  MSDHitFilter::MSDHitFilter(const Parameters& conf)
    : art::EDFilter{conf}
    , _tag(conf().tag())
    , _minHits(conf().minHits())
    , _nevt(0)
    , _npass(0)
  {
  }

  //-----------------------------------------------------------------------------
  // Check if the hit is accepted
  bool MSDHitFilter::goodHit(const MSDHit& hit) {
    if(!hit.hasTime()) return false; // time must be defined for the hit
    return true;
  }

  //-----------------------------------------------------------------------------
  bool MSDHitFilter::filter(art::Event& event) {

    // Count total events seen
    ++_nevt;

    // Filter flag
    bool passed = true;

    auto handle = event.getValidHandle<MSDHitCollection>(_tag);
    const auto hits = handle.product();
    int naccepted = 0;
    for(const auto& hit : *hits) if(goodHit(hit)) ++naccepted;

    passed &= naccepted >= _minHits;

    if (passed) ++_npass;

    // Return the result
    return passed;
  }

  //-----------------------------------------------------------------------------
  bool MSDHitFilter::endRun(art::Run& run) {
    // Print a summary of the filter results
    const float rate = (_nevt > 0) ? float(_npass)/float(_nevt) : 0.f;
    TLOG(TLVL_DEBUG + 2) << "passed " << _npass << " events out of " << _nevt << " for a ratio of " << rate;
    _nevt = 0; _npass = 0; // reset for the next run
    return true;
  }
}

DEFINE_ART_MODULE(mu2e::MSDHitFilter)
