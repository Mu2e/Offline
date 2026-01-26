#include "Offline/DataProducts/inc/FilterFraction.hh"
#include <stdexcept>

namespace mu2e {
  FilterFraction& FilterFraction::operator +=(FilterFraction const& other) {
    if(other.type() != type() || (type()<= nonominal && other.nominalFraction() != nominalFraction()))
      throw std::runtime_error("nominal filter fractions conflict");
    nseen_ += other.nSeen();
    npassed_ += other.nPassed();
    return *this;
  }
  FilterFraction FilterFraction::operator + (FilterFraction const& other) const {
    auto retval = *this;
    retval += other;
    return retval;
  }
  // concatenate with an upstream filter. The values must match!
  FilterFraction FilterFraction::chain(FilterFraction const& upstream) const {
    if(upstream.nPassed() != nSeen())throw std::runtime_error("Filter counts conflict");
    return FilterFraction(chained,-1.0,upstream.nSeen(),nPassed());
  }
}
