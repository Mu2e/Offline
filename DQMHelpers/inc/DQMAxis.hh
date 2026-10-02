#ifndef DQMHelpers_inc_DQMAxis_hh
#define DQMHelpers_inc_DQMAxis_hh
// One histogram axis, as a constant.
//
// DQM histograms are compared and merged across runs and run periods, so their
// binning must not depend on the job, the input or the geometry. Every axis is
// therefore a constexpr DQMAxis in the client header, and DQMHistSet books only
// from these -- there is no overload that takes loose numbers, and no FHiCL path
// that reaches a binning. Changing one is a code change that bumps the client's
// binning version.
//
// Original Author: R. Mina

#include <string>
#include <vector>

namespace mu2e {

struct DQMAxis {
  int n{1};
  double lo{0.};
  double hi{1.};

  constexpr DQMAxis() = default;
  constexpr DQMAxis(int nBins, double low, double high) :
      n(nBins), lo(low), hi(high)
  {}

  // Integer quantities: one bin per value, edges at the half-integers.
  static constexpr DQMAxis Counts(int first, int last)
  {
    return DQMAxis(last - first + 1, first - 0.5, last + 0.5);
  }
  // Symmetric range of a given bin width, e.g. a dt axis.
  static constexpr DQMAxis Symmetric(double range, double binWidth)
  {
    return DQMAxis(static_cast<int>(2. * range / binWidth), -range, range);
  }

  constexpr double width() const { return (hi - lo) / n; }
  bool valid() const { return n > 0 && hi > lo; }
  std::string describe() const;
};

// Variable-width bins, for the rare axis that needs them.
struct DQMVarAxis {
  std::vector<double> edges;
  bool valid() const { return edges.size() > 1; }
};

} // namespace mu2e

#endif /* DQMHelpers_inc_DQMAxis_hh */
