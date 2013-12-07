#ifndef TrackerMC_StrawElectronics_hh
#define TrackerMC_StrawElectronics_hh
//
// StrawElectronics collects the electronics response behavior of a Mu2e straw in
// several functions and parameters
//
// $Id: StrawElectronics.hh,v 1.1 2013/12/07 19:51:42 brownd Exp $
// $Author: brownd $
// $Date: 2013/12/07 19:51:42 $
//
// Original author David Brown, LBNL
//

// C++ includes
#include <iostream>
#include <array>
#include <vector>
#include <utility>

// Mu2e includes
#include "DataProducts/inc/StrawIndex.hh"
#include "RecoDataProducts/inc/StrawEnd.hh"
#include "RecoDataProducts/inc/StrawDigi.hh"
#include "TrackerMC/inc/StrawHitlet.hh"
#include "fhiclcpp/ParameterSet.h"
#include <array>
#include <vector>

namespace mu2e {
  class StrawElectronics {
    public:
      // construct from parameters
      StrawElectronics(fhicl::ParameterSet const& pset);
      virtual ~StrawElectronics();
      // analog response 
      virtual double amplifierResponse(double time,StrawHitlet const& hitlet) const; // mvolts per pCoulomb
  // digization
      virtual unsigned short adcResponse(double mvolts) const; // ADC response to analog inputs
      virtual unsigned long tdcResponse(double time) const; // TDC response to a given time
      void digitizeWaveform(std::vector<double> const& wf,StrawDigi::ADCWaveform& adc) const;
      void digitizeTimes(std::array<double,2> const& times,StrawDigi::TDCValues& tdc) const;
      bool combineEnds(double t1, double t2) const; // are times from 2 ends combined into a single digi?
// accessors
      double adcLSB() const { return _ADCLSB; } //LSB in mvolts
      double tdcLSB() const { return _TDCLSB; } //LSB in nseconds
      unsigned short maxADC() const { return _maxADC; }
      unsigned long maxTDC() const { return _maxTDC; }
      size_t nADCSamples() const { return _nADC; }
      double adcPeriod() const { return _ADCPeriod; } // period of ADC clock in nsec
      double adcOffset() const { return _ADCOffset; }
      unsigned maxDTDC() const { return _maxDTDC; } // maximum TDC difference between ends
      void adcTimes(double time, std::vector<double>& adctimes) const; // sampling times of ADC
    private:
      // scale factor between charge and voltage (milliVolts from picoCoulombs)
      double _dVdQ;
      double _trise, _tfall, _norm; // rise and fall times.  These should depend on hitlet type, FIXME!! 
      // add some noise parameter: FIXME!!!
      double _ADCLSB; // least-significant bit of ADC (mVolts)
      unsigned short _maxADC; // maximum ADC value
      size_t _nADC; // Number of ADC samples
      double _ADCPeriod; // ADC period in nsec
      double _ADCOffset; // Offset of 1st ADC sample WRT threshold crossing (nsec)
      double _TDCLSB; // least-significant bit of TDC (nsecs)
      unsigned long _maxTDC; // maximum TDC value
      unsigned _maxDTDC; // maximum TDC difference for digitizer to combine to hits in different ends to the same straw digi
  };
}

#endif

