#ifndef RecoDataProducts_STMDigi_hh
#define RecoDataProducts_STMDigi_hh
//
// STMDigi is the offline representation of the raw data readout from a
// straw.  Each STMDigi represents a single digitization, where the input
// signal to either end of the straw went over threshold.  It records the
// TDC value from both ends, and the ADC waveform train.  Parameters describing
// the conversion from digital to physical units and the waveform frequency and
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
