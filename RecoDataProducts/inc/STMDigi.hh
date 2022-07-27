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

#include "Offline/DataProducts/inc/STMTypes.hh"

namespace mu2e {

  class STMDigi {
  public:
    STMDigi() : _trigNum(0), _trigType(STMTrigType(0, 0, 0)), _trigTime(0), _trigTimeOffset(0), _extra(0), _flag(STMDigiFlag::kUnknown), _adcs(std::vector<int16_t>()){};

    STMDigi(uint32_t trigNum, STMTrigType trigType, uint64_t trigTime, uint32_t trigTimeOffset, uint32_t extra, STMDigiFlag flag, std::vector<int16_t> adcs) : _trigNum(trigNum), _trigType(trigType), _trigTime(trigTime), _trigTimeOffset(trigTimeOffset), _extra(extra), _flag(flag), _adcs(adcs) {};

    // Simpler constructor for the simulation
    STMDigi(int channel, int tdc, std::vector<int16_t> adcs) : _trigType(STMTrigType(0, channel, 0)), _trigTime(tdc), _adcs(adcs) { }

    uint32_t trigNum() const { return _trigNum; }

    STMTrigType trigType() const { return _trigType; }
    STMTriggerMode mode() const { return _trigType.mode(); }
    STMChannel channel() const { return _trigType.channel(); }
    STMDataType type() const { return _trigType.type(); }


    uint64_t trigTime() const { return _trigTime; }
    uint32_t trigTimeOffset() const { return _trigTimeOffset; }

    // information stored in extra is type-dependent
    uint16_t baselineMean() const {
      if (type() == STMDataType::kMWD) {
        return (_extra & 0xFFFF);
      }
      else {
        throw cet::exception("STMDataError") << "Trying to call baselineMean() on a digi that is not STMDataType::kMWD" << std::endl;
      }
    }

    uint16_t baselineRMS() const {
      if (type() == STMDataType::kMWD) {
        return (_extra & 0xFFFF0000);
      }
      else {
        throw cet::exception("STMDataError") << "Trying to call baselineRMS() on a digi that is not STMDataType::kMWD" << std::endl;
      }
    }


    STMDigiFlag flag() const { return _flag; }
    const std::vector<int16_t>& adcs() const { return _adcs; }

  private:
    uint32_t _trigNum;  // trigger number
    STMTrigType _trigType; // combination of: trigger mode (external (beam) or internal (source)),
                           // channel (HPGe or LaBr), and
                           // data type (unsuppressed, zero-suppressed, MWD etc.
    uint64_t _trigTime; // trigger time [ct]
    uint32_t _trigTimeOffset; // time offset from trigger to first ADC value [ct]
    uint32_t _extra;  // type-dependent (i.e. for MWD digis will be mean and RMS of baseline,for PQ digis would be quality)
    STMDigiFlag _flag;    // various error flags
    std::vector<int16_t> _adcs; // vector of ADC values for the waveform
  };

  typedef std::vector<mu2e::STMDigi> STMDigiCollection;
}
#endif
