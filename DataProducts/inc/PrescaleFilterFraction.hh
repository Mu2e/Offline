//
// art product to record the prescale filter fraction.
// Original author: Dave Brown (LBNL) 2026
//
#ifndef DataProducts_PrescaleFilterFraction_hh
#define DataProducts_PrescaleFilterFraction_hh
#include "Offline/DataProducts/inc/FilterFraction.hh"
namespace mu2e {
  class PrescaleFilterFraction : public FilterFraction {
    public:
      PrescaleFilterFraction(uint32_t prescale, uint64_t nseen, uint64_t npassed) : FilterFraction(nseen, npassed), prescale_(prescale) {}
      PrescaleFilterFraction(uint32_t prescale) : prescale_(prescale) {}
      PrescaleFilterFraction(){}
      // accessors
      uint32_t prescale() const { return prescale_; }
      double prescaleFraction() const { return 1.0/double(prescale_); }
      // concatenate multiple subruns in the same path
      PrescaleFilterFraction& operator +=(PrescaleFilterFraction const& other);
      PrescaleFilterFraction operator + (PrescaleFilterFraction const& other) const;
      // concatenate with an upstream prescale filter
      PrescaleFilterFraction chain(PrescaleFilterFraction const& upstream) const;
    private:
      uint32_t prescale_ = 0; // number of events to process for 1 to pass (on average)
  };
}
#endif
