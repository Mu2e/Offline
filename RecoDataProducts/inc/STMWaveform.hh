#ifndef RecoDataProducts_STMWaveform_hh
#define RecoDataProducts_STMWaveform_hh
//
// Data product that represents the digitized waveforms coming from the STM detectors
// This is used for both unsuppressed and zero-suppressed waveforms
//

// C++ includes
#include <iostream>
#include <vector>
#include <array>
#include <Rtypes.h>

#include "Offline/DataProducts/inc/STMChannel.hh"

namespace mu2e {

  class STMWaveform {
  public:
    STMWaveform() : _trigTimeOffset(0), _adcs(std::vector<int16_t>()){};

    STMWaveform(uint32_t trigTimeOffset, std::vector<int16_t> adcs) : _trigTimeOffset(trigTimeOffset), _adcs(adcs) {};


    uint32_t trigTimeOffset() const { return _trigTimeOffset; }
    const std::vector<int16_t>& adcs() const { return _adcs; }

  private:
    uint32_t _trigTimeOffset; // time offset from EWT? to first ADC value [ct]
    std::vector<int16_t> _adcs; // vector of ADC values for the waveform
  };

  typedef std::vector<STMWaveform> STMWaveformCollection;
}
#endif
