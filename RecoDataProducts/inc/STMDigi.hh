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
    STMDigi() {};
    STMDigi(uint16_t trigNum, uint64_t trigTime, uint64_t trigTimeOffset, uint16_t baselineMean, uint16_t baselineRMS, uint16_t nDrop, std::vector<int16_t> adcs) : _trigNum(trigNum), _trigTime(trigTime), _trigTimeOffset(trigTimeOffset), _baselineMean(baselineMean), _baselineRMS(baselineRMS), _nDrop(nDrop), _adcs(adcs) {};

    uint16_t trigNum() const { return _trigNum; }
    uint64_t trigTime() const { return _trigTime; }
    uint64_t trigTimeOffset() const { return _trigTimeOffset; }
    uint16_t baselineMean() const { return _baselineMean; }
    uint16_t baselineRMS() const { return _baselineRMS; }
    uint16_t nDrop() const { return _nDrop; }
    const std::vector<int16_t>& adcs() const { return _adcs; }

  private:
    uint16_t _trigNum;
    uint64_t _trigTime;
    uint64_t _trigTimeOffset;
    uint16_t _baselineMean;
    uint16_t _baselineRMS;
    uint16_t _nDrop;
    std::vector<int16_t> _adcs;
  };
  typedef std::vector<mu2e::STMDigi> STMDigiCollection;
}
#endif
