#include "Offline/DataProducts/inc/PrescaleFilterFraction.hh"
#include <stdexcept>

namespace mu2e {
  PrescaleFilterFraction& PrescaleFilterFraction::operator +=(PrescaleFilterFraction const& other) {
    if(other.prescale() != prescale())throw std::runtime_error("Prescale values conflict");
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
    // the net prescale of a chain is the product of the individual prescale values
    auto retval = PrescaleFilterFraction(prescale()*upstream.prescale(),upstream.nSeen(),nPassed());
    retval.setChained();
    return retval;
  }
}
