#ifndef RecoDataProducts_STMHit_hh
#define RecoDataProducts_STMHit_hh
//
// Data product represnting a time and energy of a hit in the STM
//

// C++ includes
#include <iostream>
#include <vector>
#include <array>
#include <Rtypes.h>

namespace mu2e {

  class STMHit {
  public:
    STMHit() {};
    STMHit(float time, float energy) : _time(time), _energy(energy) {};

    float time() const { return _time; }
    float energy() const { return _energy; }

  private:
    float _time;
    float _energy;
  };
  typedef std::vector<mu2e::STMHit> STMHitCollection;
}
#endif
