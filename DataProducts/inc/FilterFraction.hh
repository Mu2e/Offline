//
// art product to record the selection fraction of a filter. This product could be
// part of the raw data output (trigger prescale) or simulation presampling.
// Original author: Dave Brown (LBNL) 2026
//
#ifndef DataProducts_FilterFraction_hh
#define DataProducts_FilterFraction_hh
#include <cstdint>
namespace mu2e {
  class FilterFraction {
    public:
      FilterFraction(uint64_t nseen, uint64_t npassed) : nseen_(nseen), npassed_(npassed) {}
      // default
      FilterFraction(){}
      virtual ~FilterFraction(){}
      // accessors
      double filterFraction() const { return nseen_ > 0 ? double(npassed_)/double(nseen_) : 0.0; }
      uint64_t nSeen() const { return nseen_; }
      uint64_t nPassed() const { return npassed_; }
      bool chained() const { return chained_; }
      // concatenate multiple subruns in the same path
      FilterFraction& operator +=(FilterFraction const& other);
      FilterFraction operator + (FilterFraction const& other) const;
      // concatenate with an upstream filter
      FilterFraction chain(FilterFraction const& upstream) const;
    private:
      uint64_t nseen_ = 0; // number of events processed by this filter
      uint64_t npassed_ = 0; // number of events passed by this filter
      bool chained_ = false; // is this product the result of a chain of filters?
    protected:
      void setChained() { chained_ = true; }
  };
}
#endif
