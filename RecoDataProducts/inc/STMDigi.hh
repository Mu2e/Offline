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
    STMDigi(int tdc, int adc) : _tdc(tdc), _adc(adc) {};

    int tdc() const { return _tdc; }
    int adc() const { return _adc; }

  private:
    int _tdc;
    int _adc;
  };
  typedef std::vector<mu2e::STMDigi> STMDigiCollection;
}
#endif
