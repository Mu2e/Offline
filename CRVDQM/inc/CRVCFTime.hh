#ifndef CRVDQM_inc_CRVCFTime_hh
#define CRVDQM_inc_CRVCFTime_hh
// Constant-fraction timing for CRV waveforms. Returns the time within the
// waveform in ns; the caller adds startTDC * CRVDigitizationPeriod.
// Ported from otsdaq-mu2e-crv ArtModules/CrvCFTime.hh (mu2e/ots_ops).
//
// The threshold is a fraction of the waveform's global maximum, crossed on the
// way up to that maximum. For a waveform holding two pulses it times the
// larger one, which need not be the first.

#include "Offline/CRVConditions/inc/CRVDigitizationPeriod.hh"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace mu2e {

struct CFResult {
  double time_ns{
      std::numeric_limits<double>::quiet_NaN()}; // time within waveform [ns]
  bool valid{false};
  int16_t baseline{0};
  int16_t peak{0};
};

inline CFResult cfTime(const std::vector<int16_t>& adcs,
                       double fraction,
                       int minAmplitude,
                       double digitizationPeriod = CRVDigitizationPeriod)
{
  CFResult r;
  if (adcs.size() < 3) {
    return r;
  }

  r.baseline = adcs[0];

  auto it = std::max_element(adcs.begin(), adcs.end());
  r.peak = *it;

  double amplitude = static_cast<double>(r.peak) - r.baseline;
  if (amplitude <= 0 || amplitude < minAmplitude) {
    return r;
  }

  double threshold = r.baseline + fraction * amplitude;

  std::size_t peakIdx = static_cast<std::size_t>(std::distance(adcs.begin(), it));

  for (std::size_t i = 1; i <= peakIdx; ++i) {
    if (adcs[i] >= threshold && adcs[i - 1] < threshold) {
      double denom = static_cast<double>(adcs[i]) - adcs[i - 1];
      double frac = (denom != 0.0) ? (threshold - adcs[i - 1]) / denom : 0.0;
      r.time_ns = ((i - 1) + frac) * digitizationPeriod;
      r.valid = true;
      return r;
    }
  }

  return r;
}

} // namespace mu2e

#endif /* CRVDQM_inc_CRVCFTime_hh */
