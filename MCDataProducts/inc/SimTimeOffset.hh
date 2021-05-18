#ifndef MCDataProducts_SimTimeOffset_hh
#define MCDataProducts_SimTimeOffset_hh
//
//  Simple class to represent a global time offset, used in resampling mixing
//  Original author: David Brown (LBNL) 4/15/2021
//
#include <vector>
namespace mu2e {
  struct SimTimeOffset {
    SimTimeOffset(double toff=0.0) : timeOffset_(toff) {}
    double timeOffset_; // global event time offset
  };
  typedef std::vector<SimTimeOffset> SimTimeOffsetCollection;
}
#endif
