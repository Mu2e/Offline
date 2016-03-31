#include "GeneralUtilities/inc/NUBinning.hh"

#include <iterator>

namespace mu2e {
  const NUBinning::IndexType NUBinning::nobin(-1);

  NUBinning::IndexType NUBinning::findBin(double x) const {
    IndexType ibin = nobin;

    auto iter = std::upper_bound(binBoundaries_.begin(), binBoundaries_.end(), x);
    if((iter != binBoundaries_.begin())&&(iter != binBoundaries_.end())) {
      ibin = std::distance(binBoundaries_.begin(), iter) - 1;
    }

    return ibin;
  }

  std::ostream& operator<<(std::ostream& os, const NUBinning& b) {
    os<<"NUBinning(";
    std::copy(b.binBoundaries().begin(), b.binBoundaries().end(),  std::ostream_iterator<double>(os, " "));
    os<<")";
    return os;
  }

}
