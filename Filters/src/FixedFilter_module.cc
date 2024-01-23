//
// Trivial filter with a return value fixed by config.  This is useful for
// creating placeholders in sequences or swapping out real filters with
// pass-throughs
// original author: David Brown LBNL
//

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace mu2e {
  class FixedFilter : public art::EDFilter {
    public:
      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<bool> returnValue { Name("ReturnValue"), Comment("Fixed value returned by this filter every event") };
      };

      using Parameters = art::EDFilter::Table<Config>;
      explicit FixedFilter(const Parameters& conf);
      virtual bool filter(art::Event& event) override;
    private:
      bool retval_; // return value
  };

  FixedFilter::FixedFilter(const Parameters& conf)
    : art::EDFilter{conf}
  , retval_(conf().returnValue())
  {}

  bool FixedFilter::filter(art::Event&) {
    return retval_;
  }

}

DEFINE_ART_MODULE(mu2e::FixedFilter)
