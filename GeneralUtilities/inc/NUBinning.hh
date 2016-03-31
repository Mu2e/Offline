// A non-uniform binning representation.
//
// Andrei Gaponenko, 2016

#ifndef GeneralUtilities_inc_NUBinning_hh
#define GeneralUtilities_inc_NUBinning_hh

#include <vector>
#include <stdexcept>
#include <algorithm>
#include <ostream>

namespace mu2e {

  class NUBinning{
  public:
    typedef unsigned long IndexType;
    static const IndexType nobin; // like string::npos

    // returns nobin for out of range
    IndexType findBin(double x) const;

    IndexType nbins() const {return binBoundaries_.size() - 1;}
    const std::vector<double>& binBoundaries() const { return binBoundaries_; }

    NUBinning() = default;

    // Bin boundaries input.  The range must contain
    // at least two values (to define one bin)
    template<class Iter> NUBinning(Iter begin, Iter end)
      : binBoundaries_(begin, end)
    {
      if(binBoundaries_.size() < 2) {
        throw std::runtime_error("NUBinning(): Error: too few input to the contructor\n");
      }
      if(!std::is_sorted(binBoundaries_.begin(), binBoundaries_.end())) {
        throw std::runtime_error("NUBinning(): Error: inputs must be sorted\n");
      }
    }

  private:
    std::vector<double> binBoundaries_;
  };

  std::ostream& operator<<(std::ostream&, const NUBinning& b);
}

#endif/* GeneralUtilities_inc_NUBinning_hh */
