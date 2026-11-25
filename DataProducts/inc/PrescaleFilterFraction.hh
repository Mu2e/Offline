//
// art product to record the prescale filter fraction.
// Original author: Dave Brown (LBNL) 2026
//
// A prescale is either discrete (pass 1 in N events, N = prescale()) or fractional
// (pass each event with probability prescaleFraction(), e.g. 0.3). Discrete products are made
// with the constructors and fractional ones with fromFraction(). A disabled prescale is
// prescale() == 0 (discrete) or prescaleFraction() == 0 (fractional), and passes no events.
//
// Schema note: fraction_ was added after the first version. It defaults to a negative value,
// meaning "discrete, derive the fraction from prescale_", so products read from older files
// keep their discrete meaning without a ROOT read rule.
//
#ifndef DataProducts_PrescaleFilterFraction_hh
#define DataProducts_PrescaleFilterFraction_hh
#include "Offline/DataProducts/inc/FilterFraction.hh"
namespace mu2e {
  class PrescaleFilterFraction : public FilterFraction {
    public:
      // discrete prescale: pass 1 in prescale events
      PrescaleFilterFraction(uint32_t prescale, uint64_t nseen, uint64_t npassed) : FilterFraction(nseen, npassed), prescale_(prescale) {}
      PrescaleFilterFraction(uint32_t prescale) : prescale_(prescale) {}
      PrescaleFilterFraction(){}
      // fractional prescale: pass each event with probability fraction, in [0, 1]
      static PrescaleFilterFraction fromFraction(double fraction, uint64_t nseen = 0, uint64_t npassed = 0);
      // accessors
      bool isFractional() const { return fraction_ >= 0.0; }
      // discrete prescale factor; 0 for a fractional prescale (or a disabled discrete one)
      uint32_t prescale() const { return isFractional() ? 0 : prescale_; }
      // expected fraction of events passing the prescale, for either kind; 0 if disabled
      double prescaleFraction() const;
      // concatenate multiple subruns in the same path; throws unless both are the same kind with the same prescale
      PrescaleFilterFraction& operator +=(PrescaleFilterFraction const& other);
      PrescaleFilterFraction operator + (PrescaleFilterFraction const& other) const;
      // concatenate with an upstream prescale filter
      PrescaleFilterFraction chain(PrescaleFilterFraction const& upstream) const;
    private:
      uint32_t prescale_ = 0; // number of events to process for 1 to pass (on average); discrete only
      double fraction_ = -1.0; // pass probability for a fractional prescale; < 0 means discrete
  };
}
#endif
