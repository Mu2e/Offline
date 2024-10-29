#ifndef RecoDataProducts_CrvRecoPulseFlags_hh
#define RecoDataProducts_CrvRecoPulseFlags_hh

#include <bitset>

namespace mu2e
{
  enum CrvRecoPulseFlagEnums{failedFit=0, duplicateNoFitPulse=1, separatedDoublePulse=2, zeroNdf=3, noCalibConstPulseArea=4, noCalibConstPulseHeight=5};

  static constexpr int nFlags=8;

  typedef std::bitset<nFlags> CrvRecoPulseFlags;
}

#endif /* RecoDataProducts_CrvRecoPulseFlags_hh */

