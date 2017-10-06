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
#include <vector>
#include <utility>

// Mu2e includes
#include "Mu2eInterfaces/inc/ConditionsEntity.hh"
#include "TrackerConditions/inc/Types.hh"
#include "fhiclcpp/ParameterSet.h"

namespace mu2e {
  struct WireDistancePoint{
    double _distance;
    double _mean;
    double _normalization;
    double _sigma;
    double _t0;
    double _tmax[TrkTypes::npaths]; // time at which value is maximum
    double _linmax[TrkTypes::npaths]; // linear response to unit charge at maximum
    std::vector<double> _currentPulse;
    std::vector<double> _preampResponse;
    std::vector<double> _adcResponse;
    std::vector<double> _preampToAdc1Response;
    WireDistancePoint(){};
    WireDistancePoint(double distance, double mean, double normalization, double sigma, double t0) : 
      _distance(distance), _mean(mean), _normalization(normalization), _sigma(sigma), _t0(t0) {};
  };

  class StrawElectronics : virtual public ConditionsEntity {
    public:
      // construct from parameters
      StrawElectronics(fhicl::ParameterSet const& pset);
      virtual ~StrawElectronics();
      // linear response to a charge pulse.  This does NOT include saturation effects,
      // since those are cumulative and cannot be computed for individual charges
      double linearResponse(TrkTypes::Path ipath, double time,double charge,double distance,bool forsaturation=false) const; // mvolts per pCoulomb
      double adcImpulseResponse(double time, double charge) const;
      // Given a (linear) total voltage, compute the saturated voltage
      double saturatedResponse(double lineearresponse) const;
      // relative time when linear response is maximal
      double maxResponseTime(TrkTypes::Path ipath,double distance) const;
  // digization
      unsigned short adcResponse(double mvolts) const; // ADC response to analog inputs
      unsigned long tdcResponse(double time) const; // TDC response to a given time
      void digitizeWaveform(TrkTypes::ADCVoltages const& wf,TrkTypes::ADCWaveform& adc) const; // digitize an array of voltages at the ADC
      void digitizeTimes(TrkTypes::TDCTimes const& times,TrkTypes::TDCValues& tdc) const;
      bool combineEnds(double t1, double t2) const; // are times from 2 ends combined into a single digi?
  // interpretation of digital data
      void tdcTimes(TrkTypes::TDCValues const& tdc, TrkTypes::TDCTimes& times) const;
      double adcVoltage(unsigned short adcval) const; // mVolts
      double adcCurrent(unsigned short adcval) const; // microAmps
// accessors
      double adcLSB() const { return _ADCLSB; } //LSB in mvolts
      double tdcLSB() const { return _TDCLSB; } //LSB in nseconds
      double totLSB() const { return _TOTLSB; } //LSB in nseconds
      unsigned short maxADC() const { return _maxADC; }
      unsigned long maxTDC() const { return _maxTDC; }
      unsigned short maxTOT() const { return _maxTOT; }
      unsigned short ADCPedestal() const { return _ADCped; };
      size_t nADCSamples() const { return _nADC; }
      size_t nADCPreSamples() const { return _nADCpre; }
      double adcPeriod() const { return _ADCPeriod; } // period of ADC clock in nsec
      double adcOffset() const { return _ADCOffset; } // offset WRT clock edge for digitization
      double flashStart() const { return _flashStart; } // time flash blanking starts
      double flashEnd() const { return _flashEnd; } // time flash blanking ends
      void adcTimes(double time, TrkTypes::ADCTimes& adctimes) const; // given crossing time, fill sampling times of ADC CHECK THIS IS CORRECT IN DRAC FIXME!
      double saturationVoltage() const { return _vsat; }
      double threshold() const { return _vthresh; }
      double analogNoise(TrkTypes::Path ipath) const { return _analognoise[ipath]; }  // incoherent noise
      double strawNoise() const { return _snoise;} // coherent part of threshold circuit noise
      double deadTimeAnalog() const { return _tdeadAnalog; }
      double deadTimeDigital() const { return _tdeadDigital; }
      double clockStart() const { return _clockStart; }
      double clockJitter() const { return _clockJitter; }
      double currentToVoltage(TrkTypes::Path ipath) const { return _dVdI[ipath]; }
      double maxLinearResponse(TrkTypes::Path ipath,double distance,double charge=1.0) const;
      double peakMinusPedestalEnergyScale() const { return _pmpEnergyScale; }
      double normalization(TrkTypes::Path ipath) const { return 1.;} //FIXME
      double fallTime(TrkTypes::Path ipath) const { return 22.;} //FIXME

      void calculateResponse(std::vector<double> &poles, std::vector<double> &zeros, std::vector<double> &input, std::vector<double> &response, double dVdI);
      double truncationTime(TrkTypes::Path ipath) const { return _ttrunc[ipath];}
      double saturationTimeStep() const { return _saturationSampleFactor/_sampleRate;}
      static double _pC_per_uA_ns; // unit conversion from pC/ns to microAmp
    private:
    // generic waveform parameters
      double _dVdI[TrkTypes::npaths]; // scale factor between charge and voltage (milliVolts/picoCoulombs)
      double _ttrunc[TrkTypes::npaths]; // time to truncate signal to 0
// threshold path parameters
      double _tdeadAnalog; // electronics dead time
      double _tdeadDigital; // electronics readout dead time
      // scale factor between current and voltage (milliVolts per microAmps)
      double _vsat; // saturation parameters.  _vmax is maximum output, _vsat is where saturation starts
      double _vthresh; // threshold voltage for electronics discriminator (mVolt)
      double _snoise; // straw noise at threshold
      double _analognoise[TrkTypes::npaths]; //noise (mVolt) from the straw itself
      double _ADCLSB; // least-significant bit of ADC (mVolts)
      unsigned short _maxADC; // maximum ADC value
      unsigned short _ADCped; // ADC pedestal (reading for 0 volts)
      size_t _nADC,_nADCpre; // Number of ADC samples, presamples
      double _ADCPeriod; // ADC period in nsec
      double _ADCOffset; // Offset of 1st ADC sample WRT threshold crossing (nsec)
      double _TDCLSB; // least-significant bit of TDC (nsecs)
      unsigned long _maxTDC; // maximum TDC value
      double _TOTLSB; // least-significant bit of TOT (nsecs)
      unsigned short _maxTOT; // maximum TOT value
      double _clockStart, _clockJitter; // time TDC clock starts, and its error (common to both ends!!)
      double _flashStart, _flashEnd; // flash blanking period (no digitizations during this time!!!)
      double _pmpEnergyScale; // fudge factor for peak minus pedestal energy method
  // helper functions
      static inline double mypow(double,unsigned);

      int _responseBins;
      double _sampleRate;
      int _saturationSampleFactor;
      std::vector<double> _preampPoles;
      std::vector<double> _preampZeros;
      std::vector<double> _adcPoles;
      std::vector<double> _adcZeros;
      std::vector<double> _preampToAdc1Poles;
      std::vector<double> _preampToAdc1Zeros;
      std::vector<double> _preampToAdc2Poles;
      std::vector<double> _preampToAdc2Zeros;

      std::vector<double> _wireDistances;
      std::vector<double> _currentMeans;
      std::vector<double> _currentNormalizations;
      std::vector<double> _currentSigmas;
      std::vector<double> _currentT0s;

      std::vector<double> _currentImpulse;
      std::vector<double> _preampToAdc2Response;
      
      std::vector<WireDistancePoint> _wPoints;

  };
}

#endif

