#ifndef RecoDataProducts_STMWaveformDigi_hh
#define RecoDataProducts_STMWaveformDigi_hh
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

  class STMWaveformDigi {
  public:

    STMWaveformDigi() : _DetID(-1), _trigTimeOffset(0), _adcs(std::vector<int16_t>()), _peakpos(0){};

    STMWaveformDigi(int DetID, uint32_t trigTimeOffset, std::vector<int16_t> adcs, int peakpos) : _DetID(DetID), _trigTimeOffset(trigTimeOffset), _adcs(adcs), _peakpos(peakpos) {};
    STMWaveformDigi(int DetID, uint32_t trigTimeOffset, std::vector<int16_t> adcs) : _DetID(DetID), _trigTimeOffset(trigTimeOffset), _adcs(adcs), _peakpos(0) {};
    STMWaveformDigi(uint32_t trigTimeOffset, std::vector<int16_t> adcs) : _trigTimeOffset(trigTimeOffset), _adcs(adcs) {};

    int                     DetID()   const {return _DetID;}
    uint32_t trigTimeOffset() const { return _trigTimeOffset; }
    const std::vector<int16_t>& adcs() const { return _adcs; }
    int                     peakpos()  const {return _peakpos;}

  private:
    int _DetID;
    uint32_t _trigTimeOffset; // time offset from EWT? to first ADC value [ct]
    std::vector<int16_t> _adcs; // vector of ADC values for the waveform
    int _peakpos;
  };

  typedef std::vector<STMWaveformDigi> STMWaveformDigiCollection;
}
#endif
