#ifndef RecoDataProducts_CrvRecoPulseFlags_hh
#define RecoDataProducts_CrvRecoPulseFlags_hh

#include <bitset>

namespace mu2e 
{
  enum CrvRecoPulseFlagEnums{failedFit=0, unused=1, unused2=2};

  typedef std::bitset<8> CrvRecoPulseFlags;
}

#endif /* RecoDataProducts_CrvRecoPulseFlags_hh */

