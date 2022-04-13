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
/*
  struct Test {
    uint_16t trigNum;
 
  }
*/
  
  class STMDigi {
  public:
    STMDigi() {};
    STMDigi(uint16_t trigNum, uint32_t trigTime, uint16_t trigTimeOffset, uint16_t baselineMean, uint16_t baselineRMS, uint16_t nDrop, int adc0) : _trigNum(trigNum), _trigTime(trigTime), _trigTimeOffset(trigTimeOffset), _baselineMean(baselineMean), _baselineRMS(baselineRMS), _nDrop(nDrop), _adc0(adc0) {};
    STMDigi(uint32_t trigTime, int adc0): _trigTime(trigTime), _adc0(adc0) {};

    uint16_t trigNum() const { return _trigNum; }
    uint32_t trigTime() const { return _trigTime; }
    uint16_t trigTimeOffset() const { return _trigTimeOffset; }
    uint16_t baselineMean() const { return _baselineMean; }
    uint16_t baselineRMS() const { return _baselineRMS; }
    uint16_t nDrop() const { return _nDrop; }
    int adc0() const { return _adc0; }

  private:
    uint16_t _trigNum;
    uint32_t _trigTime;
    uint16_t _trigTimeOffset;
    uint16_t _baselineMean;
    uint16_t _baselineRMS;
    uint16_t _nDrop;
    int _adc0; 
  };
  typedef std::vector<mu2e::STMDigi> STMDigiCollection;
}
#endif
