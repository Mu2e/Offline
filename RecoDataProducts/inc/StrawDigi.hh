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
// $Id: StrawDigi.hh,v 1.3 2013/12/12 19:07:56 brownd Exp $
// $Author: brownd $
// $Date: 2013/12/12 19:07:56 $
//
// Original author David Brown, LBNL
//

// C++ includes
#include <iostream>
#include <vector>
#include <array>

// Mu2e includes
#include "DataProducts/inc/StrawIndex.hh"
#include "TrackerConditions/inc/StrawEnd.hh"
#include "TrackerConditions/inc/Types.hh"
// unfortunately the persistence requires this array dimension be
// instantiated by hand in the classesdef, so this this #define
// can't be changed without changing that too
namespace mu2e {
  
  class StrawDigi {
    public:
      StrawDigi();
      StrawDigi( StrawIndex sid, TrkTypes::TDCValues tdc, TrkTypes::TOTValues tot, TrkTypes::ADCWaveform const& adc);
      StrawDigi(StrawDigi const& other);
      StrawDigi& operator=(StrawDigi const& other);
      StrawIndex strawIndex() const { return _strawIndex; }
// TDC data are indexed according to straw end
      unsigned long TDC(StrawEnd end) const { return _tdc[end]; }
      unsigned long TOT(StrawEnd end) const { return _tot[end]; }
      TrkTypes::TDCValues const& TDC() const { return _tdc; }
      TrkTypes::TOTValues const& TOT() const { return _tot; }
      TrkTypes::ADCWaveform const& adcWaveform() const { return _adc; }
    private:
      StrawIndex  _strawIndex;      // Straw index
      TrkTypes::TDCValues _tdc; // TDC values for each end
      TrkTypes::TOTValues _tot;  // TOT values for each end
      TrkTypes::ADCWaveform _adc; // ADC waveform (sum of both ends)
  };
  typedef std::vector<mu2e::StrawDigi> StrawDigiCollection;
}
#endif
