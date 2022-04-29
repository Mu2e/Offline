#ifndef RecoDataProducts_STMDigi_hh
#define RecoDataProducts_STMDigi_hh
//
// Data product that represents the digitized signal coming from the STM
//

// C++ includes
#include <iostream>
#include <vector>
#include <array>
#include <Rtypes.h>

namespace mu2e {

  class STMDigi {
  public:
    STMDigi() : _trigNum(0), _trigMode(0), _channel(0), _trigTime(0), _trigTimeOffset(0), _baselineMean(0), _baselineRMS(0), _nDrop(0), _adcs(std::vector<int16_t>()){};

    STMDigi(uint32_t trigNum, uint16_t trigMode, uint16_t channel, uint64_t trigTime, uint32_t trigTimeOffset, uint16_t baselineMean, uint16_t baselineRMS, uint16_t nDrop, std::vector<int16_t> adcs) : _trigNum(trigNum), _trigMode(trigMode), _channel(channel), _trigTime(trigTime), _trigTimeOffset(trigTimeOffset), _baselineMean(baselineMean), _baselineRMS(baselineRMS), _nDrop(nDrop), _adcs(adcs) {};

    // Simpler constructor for the simulation
    STMDigi(int tdc, int adc) : _trigTime(tdc), _adcs(std::vector<int16_t>()) {
      _adcs.push_back(adc);
    }

    uint32_t trigNum() const { return _trigNum; }
    uint16_t trigMode() const { return _trigMode; }
    uint16_t channel() const { return _channel; }
    uint64_t trigTime() const { return _trigTime; }
    uint32_t trigTimeOffset() const { return _trigTimeOffset; }
    uint16_t baselineMean() const { return _baselineMean; }
    uint16_t baselineRMS() const { return _baselineRMS; }
    uint16_t nDrop() const { return _nDrop; }
    const std::vector<int16_t>& adcs() const { return _adcs; }

  private:
    uint32_t _trigNum;  // trigger number
    uint16_t _trigMode; // external (beam) or internal (source)
    uint16_t _channel;  // HPGe or LaBr
    uint64_t _trigTime; // trigger time [ct]
    uint32_t _trigTimeOffset; // time offset from trigger to first ADC value [ct]
    uint16_t _baselineMean;   // mean baseline (calculated by MWD algorithm)
    uint16_t _baselineRMS;    // RMS of baseline (calculated by MWD algorithm)
    uint16_t _nDrop;    // number of dropped packets
    std::vector<int16_t> _adcs; // vector of ADC values for the waveform
  };

  typedef std::vector<mu2e::STMDigi> STMDigiCollection;
}
#endif
