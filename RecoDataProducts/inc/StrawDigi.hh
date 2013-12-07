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
// $Id: StrawDigi.hh,v 1.1 2013/12/07 19:50:42 brownd Exp $
// $Author: brownd $
// $Date: 2013/12/07 19:50:42 $
//
// Original author David Brown, LBNL
//

// C++ includes
#include <iostream>
#include <vector>
//#include <array>

// Mu2e includes
#include "DataProducts/inc/StrawIndex.hh"
#include "RecoDataProducts/inc/StrawEnd.hh"
// unfortunately the persistence requires this array dimension be
// instantiated by hand in the classesdef, so this this #define
// can't be changed without changing that too
//#define NADC 10
namespace mu2e {
  
  class StrawDigi {
    public:
// delimit the range of physical channels
      enum TDCChannel {zero=0,one=1};
// unfortunately genreflex doesn't understand array, so I can't use these for now
//      typedef std::array<unsigned short,NADC> ADCWaveform;
//      typedef std::array<unsigned long,2> TDCValues;
      typedef std::vector<unsigned short> ADCWaveform;
      typedef unsigned long TDCValues[2];
      StrawDigi();
      StrawDigi( StrawIndex sid, TDCValues tdc, ADCWaveform const& adc);
      StrawDigi(StrawDigi const& other);
      StrawDigi& operator=(StrawDigi const& other);
      StrawIndex strawIndex() const { return _strawIndex; }
// TDC data are indexed according to electronics channel
      unsigned long TDC(TDCChannel channel) { return _tdc[(size_t)channel]; }
// TDC according to the physical straw end.  This requires accessing conditions data
      unsigned long TDC(StrawEnd end) const;
      ADCWaveform const& adcWaveform() const { return _adc; }
    private:
      StrawIndex  _strawIndex;      // Straw index
      TDCValues _tdc;
      ADCWaveform _adc;
  };
}
#endif
