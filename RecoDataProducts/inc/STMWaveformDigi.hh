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
    // Initialise all variables
    STMWaveformDigi() : _DetID(0), _EWT(0), _DTCtime(0), _ADCtime(0), _trigTimeOffset(0), _peak_fitTime1(0), _peak_fitTime2(0), _peak_sep(0), _adcs(std::vector<int16_t>()) {};
    // Constructor for timing plus trig offset
    STMWaveformDigi(int16_t DetID, uint64_t EWT, uint64_t DTCtime, uint64_t ADCtime, uint32_t trigTimeOffset, std::vector<int16_t> &adcs) : _DetID(DetID), _EWT(EWT), _DTCtime(DTCtime), _ADCtime(ADCtime), _trigTimeOffset(trigTimeOffset), _adcs(adcs) {};
    // Constructor for peak fitting
    STMWaveformDigi(int16_t DetID, uint64_t EWT, uint64_t DTCtime, uint64_t ADCtime, uint32_t trigTimeOffset, double peak_fitTime1, double peak_fitTime2, double peak_sep, std::vector<int16_t> &adcs) : _DetID(DetID), _EWT(EWT), _DTCtime(DTCtime), _ADCtime(ADCtime), _trigTimeOffset(trigTimeOffset), _peak_fitTime1(peak_fitTime1), _peak_fitTime2(peak_fitTime2), _peak_sep(peak_sep), _adcs(adcs) {};
    // Basic constructor
    STMWaveformDigi(uint32_t trigTimeOffset, std::vector<int16_t> &adcs) : _DetID(0), _EWT(0), _DTCtime(0), _ADCtime(0), _trigTimeOffset(trigTimeOffset), _peak_fitTime1(0), _peak_fitTime2(0), _peak_sep(0), _adcs(adcs) {};

    int16_t                     DetID  () const { return _DetID; }
    uint64_t                    EWT    () const { return _EWT; }
    uint64_t                    DTCtime() const { return _DTCtime; }
    uint64_t                    ADCtime() const { return _ADCtime; }
    uint32_t                    trigTimeOffset() const { return _trigTimeOffset; }
    double                      peak_fitTime1 () const { return _peak_fitTime1; }
    double                      peak_fitTime2 () const { return _peak_fitTime2; }
    double                      peak_sep      () const { return _peak_sep; }
    const std::vector<int16_t>& adcs   () const { return _adcs; }
    
  private:
    int16_t  _DetID;
    uint64_t _EWT;
    uint64_t _DTCtime;
    uint64_t _ADCtime;
    uint32_t _trigTimeOffset; // time offset from EWT? to first ADC value [ct]
    double   _peak_fitTime1; // fit time of first rising edge (ns)
    double   _peak_fitTime2; // fit time of second rising edge (ns)
    double   _peak_sep; // separation time (ns)
    std::vector<int16_t> _adcs; // vector of ADC values for the waveform
  };

  typedef std::vector<STMWaveformDigi> STMWaveformDigiCollection;
}
#endif
