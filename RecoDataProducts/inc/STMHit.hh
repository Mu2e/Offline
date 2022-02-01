#ifndef RecoDataProducts_STMHit_hh
#define RecoDataProducts_STMHit_hh
//
// STMHit is the offline representation of the raw data readout from a
// straw.  Each STMHit represents a single hittization, where the input
// signal to either end of the straw went over threshold.  It records the
// TDC value from both ends, and the ADC waveform train.  Parameters describing
// the conversion from hittal to physical units and the waveform frequency and
// phase are stored in conditions objects.
//
//
// Original author David Brown, LBNL
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
