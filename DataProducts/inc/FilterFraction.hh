//
// Run-level art product to record the selection fraction of a filter. This product could be
// part of the raw data output (trigger prescale) or simulation.
// Original author: Dave Brown (LBNL) 2026
//
#ifndef DataProducts_FilterFraction_hh
#define DataProducts_FilterFraction_hh
namespace mu2e {
  class FilterFraction {
    public:
      enum FilterType { constant=0,sculpted}; // sculpted means the selection depends on the event data content, not just purely (pseudo) random.
      FilterFraction(FilterType type, double nominal, double actual) : type_(type),nominal_(nominal), actual_(actual) {}
      // accessors
      FilterType type() const { return type_; }
      double nominalFraction() const { return nominal_; }
      double actualFraction() const { return actual_; }
    private:
      FilterType type_; // type of filtering performed
      double nominal_; // nominal fraction (<1) of events expected to be kept by the filter
      double actual_; // fraction (<1) of events actually kept by the filter
  };
}
#endif
