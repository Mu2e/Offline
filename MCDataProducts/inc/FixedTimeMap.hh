// A data product that allows you to chose proton times at generator
// time instead of after Geant4 simulation. All sim particles from this
// event will be given the same proton time

#ifndef FixedProtonPulseTimeAssns_hh
#define FixedProtonPulseTimeAssns_hh

#include <map>

#include "canvas/Persistency/Common/Ptr.h"


namespace mu2e {
  class FixedTimeMap {
  public:
    double time() const { return time_; }

    void SetTime(double t){ time_ = t;}

    // defautl ctr required by ROOT persistency
    FixedTimeMap() : time_(-1) {}

  private:
    double time_;
  };

}

#endif/*FixedProtonPulseTimeAssns_hh*/
