#ifndef RecoDataProducts_StrawDigi_hh
#define RecoDataProducts_StrawDigi_hh
//
// StrawDigi is the offline representation of the raw data readout from a
// straw.  Each StrawDigi represents a single digitization, where the input
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

#include "canvas/Persistency/Common/Ptr.h"

// Mu2e includes
#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/StrawEnd.hh"
#include "RecoDataProducts/inc/StrawDigiFlag.hh"
#include "DataProducts/inc/TrkTypes.hh"
// unfortunately the persistence requires this array dimension be
// instantiated by hand in the classesdef, so this this #define
// can't be changed without changing that too
namespace mu2e {

  class StrawDigiADCWaveform {
    public:
      StrawDigiADCWaveform() = default;
      StrawDigiADCWaveform(TrkTypes::ADCWaveform const& adc) : _adc(adc) {};

      TrkTypes::ADCWaveform const& samples() const { return _adc; }

    private:
      TrkTypes::ADCWaveform _adc;
  };
  typedef std::vector<mu2e::StrawDigiADCWaveform> StrawDigiADCWaveformCollection;
  
  class StrawDigi {
    public:
      StrawDigi() = default;
      StrawDigi(StrawId sid, TrkTypes::TDCValues tdc, TrkTypes::TOTValues tot, TrkTypes::ADCValue pmp);

      StrawId strawId() const { return _strawid; }
// TDC data are indexed according to straw end
      unsigned long TDC(StrawEnd end) const { return _tdc[end]; }
      unsigned long TOT(StrawEnd end) const { return _tot[end]; }
      TrkTypes::TDCValues const& TDC() const { return _tdc; }
      TrkTypes::TOTValues const& TOT() const { return _tot; }
      TrkTypes::ADCValue const& PMP() const { return _pmp; }
      StrawDigiFlag const& digiFlag() const { return _flag; }
      StrawDigiFlag& digiFlag() { return _flag; }
    private:
      StrawId  _strawid;      // Straw id
      TrkTypes::TDCValues _tdc; // TDC values for each end
      TrkTypes::TOTValues _tot;  // TOT values for each end
      StrawDigiFlag _flag; // bit flags for this digi
      TrkTypes::ADCValue _pmp; // firmware calculated peak minus pedestal
  };
  typedef std::vector<mu2e::StrawDigi> StrawDigiCollection;
}
#endif
