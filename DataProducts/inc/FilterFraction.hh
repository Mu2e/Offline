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
      enum FilterType { constant=0, nonominal, chained, unknown};
      // chained means multiple filtering steps have been chained together
      FilterFraction(FilterType type, double nominal, uint64_t nseen, uint64_t npassed) :
        type_(type),nominal_(nominal), nseen_(nseen), npassed_(npassed) {}
      // use this constructor when there is no nominal selection fraction
      FilterFraction(uint64_t nseen, uint64_t npassed) :
        type_(nonominal),nominal_(-1.0), nseen_(nseen), npassed_(npassed) {}
      // default
      FilterFraction(){}
      // accessors
      FilterType type() const { return type_; }
      bool hasNominalValue() const { return type() == nonominal; }
      double nominalFraction() const { return nominal_; }
      double actualFraction() const { return nseen_ > 0 ? double(npassed_)/double(nseen_) : 0.0; }
      uint64_t nSeen() const { return nseen_; }
      uint64_t nPassed() const { return npassed_; }
      // concatenate multiple subruns in the same path
      FilterFraction& operator +=(FilterFraction const& other);
      FilterFraction operator + (FilterFraction const& other) const;
      // concatenate with an upstream filter. The values must match!
      FilterFraction chain(FilterFraction const& upstream) const;
    private:
      FilterType type_ = unknown; // type of filtering performed
      double nominal_ = -1; // nominal fraction (<1) of events expected to be kept by the filter
      uint64_t nseen_ = 0; // number of events processed by this filter
      uint64_t npassed_ = 0; // number of events passed by this filter
  };
}
#endif
