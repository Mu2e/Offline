#ifndef STMConditions_STMEnergyCorr_hh
#define STMConditions_STMEnergyCorr_hh

// parameters used to correct digis to energy for an STM channel

namespace mu2e {

struct STMEnergyCorr {
  double p0;
  double p1;
  double p2;
};

}  // namespace mu2e

#endif
