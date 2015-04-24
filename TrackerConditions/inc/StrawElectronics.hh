#ifndef TrackerConditions_StrawElectronics_hh
#define TrackerConditions_StrawElectronics_hh
//
// StrawElectronics collects the electronics response behavior of a Mu2e straw in
// several functions and parameters
//
// $Id: StrawElectronics.hh,v 1.6 2014/03/11 16:18:01 brownd Exp $
// $Author: brownd $
// $Date: 2014/03/11 16:18:01 $
//
// Original author David Brown, LBNL
//

// C++ includes
#include <iostream>
#include <array>
#include <vector>
#include <utility>

// Mu2e includes
#include "Mu2eInterfaces/inc/ConditionsEntity.hh"
#include "fhiclcpp/ParameterSet.h"
#include <array>
#include <vector>

namespace mu2e {
  class StrawElectronics : virtual public ConditionsEntity {
    public:
// separately describe the 2 analog paths
      enum path{thresh=0,adc};
// these are copied from StrawDigi, but I don't want a direct dependency
      typedef unsigned long TDCValues[2];
      typedef std::vector<unsigned short> ADCWaveform;
 
      // construct from parameters
      StrawElectronics(fhicl::ParameterSet const& pset);
      virtual ~StrawElectronics();
      // linear response to a charge pulse.  This does NOT include saturation effects,
      // since those are cumulative and cannot be computed for individual charges
      double linearResponse(path ipath, double time,double charge) const; // mvolts per pCoulomb
      // Given a (linear) total voltage, compute the saturated voltage
      double saturatedResponse(double lineearresponse) const;
      // relative time when linear response is maximal
      double maxResponseTime(path ipath) const { return _tmax[ipath]; }
  // digization
      unsigned short adcResponse(double mvolts) const; // ADC response to analog inputs
      unsigned long tdcResponse(double time) const; // TDC response to a given time
      void digitizeWaveform(std::vector<double> const& wf,ADCWaveform& adc) const;
      void digitizeTimes(std::array<double,2> const& times,TDCValues& tdc) const;
      bool combineEnds(double t1, double t2) const; // are times from 2 ends combined into a single digi?
  // interpretation of digital data
      void tdcTimes(TDCValues const& tdc, std::array<double,2>& times) const;
      double adcVoltage(unsigned short adcval) const; // mVolts
      double adcCurrent(unsigned short adcval) const; // microAmps
// accessors
      double adcLSB() const { return _ADCLSB; } //LSB in mvolts
      double tdcLSB() const { return _TDCLSB; } //LSB in nseconds
      unsigned short maxADC() const { return _maxADC; }
      unsigned long maxTDC() const { return _maxTDC; }
      unsigned short ADCPedestal() const { return _ADCped; };
      size_t nADCSamples() const { return _nADC; }
      size_t nADCPreSamples() const { return _nADCpre; }
      double adcPeriod() const { return _ADCPeriod; } // period of ADC clock in nsec
      double adcOffset() const { return _ADCOffset; } // offset WRT clock edge for digitization
      double flashStart() const { return _flashStart; } // time flash blanking starts
      double flashEnd() const { return _flashEnd; } // time flash blanking ends
      void adcTimes(double time, std::vector<double>& adctimes) const; // sampling times of ADC
      double saturationVoltage() const { return _vsat; }
      double maximumVoltage() const { return _vmax; }
      double fallTime(path ipath) const { return _tau[ipath]; }
      double dispersion(double dlen) const { return _disp*dlen; } // dispersion width is linear in propagation length
      double threshold() const { return _vthresh; }
      double thresholdNoise() const { return _vthreshnoise; }
      double deadTime() const { return _tdead; }
      double clockStart() const { return _clockStart; }
      double clockJitter() const { return _clockJitter; }
      double currentToVoltage(path ipath) const { return _dVdI[ipath]; }
      double normalization(path ipath) const { return _norm[ipath]; }
      double maxLinearResponse(path ipath,double charge=1.0) const { return _linmax[ipath]*charge; }
      static double _pC_per_uA_ns; // unit conversion from pC/ns to microAmp
    private:
    // generic waveform parameters
      double _dVdI[2]; // scale factor between charge and voltage (milliVolts/picoCoulombs)
      double _tmax[2]; // time at which value is maximum
      double _ttrunc[2]; // time to truncate signal to 0
      double _linmax[2]; // linear response to unit charge at maximum
      double _norm[2]; // normalization factor
// adc path parameters
      double _tau[2], _freq[2]; // shaping time and associated frequency
// threshold path parameters
      double _tband, _voff, _toff;
      double _tdead; // electronics dead time
      // scale factor between current and voltage (milliVolts per microAmps)
      double _vmax, _vsat, _vdiff; // saturation parameters.  _vmax is maximum output, _vsat is where saturation starts
      double _disp; // dispersion in ns/mm;
      double _vthresh; // threshold voltage for electronics discriminator (mVolt)
      double _vthreshnoise; // threshold voltage noise width (mVolt)
// add some noise parameter: FIXME!!!
      double _ADCLSB; // least-significant bit of ADC (mVolts)
      unsigned short _maxADC; // maximum ADC value
      unsigned short _ADCped; // ADC pedestal (reading for 0 volts)
      size_t _nADC,_nADCpre; // Number of ADC samples, presamples
      double _ADCPeriod; // ADC period in nsec
      double _ADCOffset; // Offset of 1st ADC sample WRT threshold crossing (nsec)
      double _TDCLSB; // least-significant bit of TDC (nsecs)
      unsigned long _maxTDC; // maximum TDC value
      double _clockStart, _clockJitter; // time TDC clock starts, and its error (common to both ends!!)
      double _flashStart, _flashEnd; // flash blanking period (no digitizations during this time!!!)
  // helper functions
      static inline double mypow(double,unsigned);
  };
}

#endif

