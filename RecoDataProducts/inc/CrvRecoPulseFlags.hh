#ifndef RecoDataProducts_CrvRecoPulseFlags_hh
#define RecoDataProducts_CrvRecoPulseFlags_hh

#include <bitset>

namespace mu2e
{
  enum CrvRecoPulseFlagEnums{failedFit=0, duplicateNoFitPulse=1, separatedDoublePulse=2, zeroNdf=3};

  typedef std::bitset<8> CrvRecoPulseFlags;
}

#endif /* RecoDataProducts_CrvRecoPulseFlags_hh */

