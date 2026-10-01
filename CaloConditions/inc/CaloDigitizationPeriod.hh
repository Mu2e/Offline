#ifndef CaloConditions_CaloDigitizationPeriod_hh
#define CaloConditions_CaloDigitizationPeriod_hh

namespace mu2e {

  // DIRAC ADC sampling period. CaloDigi::t0() is in units of this period (digitizer clock ticks).
  constexpr double CaloDigitizationPeriod = 5.0; // ns

}  // namespace mu2e

#endif
