#include "Offline/DataProducts/inc/PrescaleFilterFraction.hh"
#include <stdexcept>
#include <string>

namespace mu2e {
  PrescaleFilterFraction PrescaleFilterFraction::fromFraction(double fraction, uint64_t nseen, uint64_t npassed) {
    if(!(fraction >= 0.0 && fraction <= 1.0))
      throw std::invalid_argument("PrescaleFilterFraction fraction must be in [0, 1], got " + std::to_string(fraction));
    PrescaleFilterFraction retval(0u, nseen, npassed);
    retval.fraction_ = fraction;
    return retval;
  }

  double PrescaleFilterFraction::prescaleFraction() const {
    if(isFractional()) return fraction_;
    return prescale_ > 0 ? 1.0/double(prescale_) : 0.0;
  }

  PrescaleFilterFraction& PrescaleFilterFraction::operator +=(PrescaleFilterFraction const& other) {
    // both products must be the same kind with the same prescale (N if discrete, fraction if fractional)
    const bool conflict = (isFractional() != other.isFractional())
      || (isFractional() ? other.fraction_ != fraction_ : other.prescale_ != prescale_);
    if(conflict) throw std::runtime_error("Prescale values conflict");
    // invoke base class
    (static_cast<FilterFraction*>(this))->operator+=(other);
    return *this;
  }
  PrescaleFilterFraction PrescaleFilterFraction::operator + (PrescaleFilterFraction const& other) const {
    auto retval = *this;
    retval += other;
    return retval;
  }
  // concatenate with an upstream filter
  PrescaleFilterFraction PrescaleFilterFraction::chain(PrescaleFilterFraction const& upstream) const {
    if(upstream.nPassed() != nSeen())throw std::runtime_error("Filter counts conflict");
    // the net prescale of a chain is the product of the individual prescales: an integer
    // product if both are discrete, otherwise the product of the pass fractions
    auto retval = (!isFractional() && !upstream.isFractional())
      ? PrescaleFilterFraction(prescale_*upstream.prescale_,upstream.nSeen(),nPassed())
      : fromFraction(prescaleFraction()*upstream.prescaleFraction(),upstream.nSeen(),nPassed());
    retval.setChained();
    return retval;
  }
}
