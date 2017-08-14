//
// StrawElectronics collects the electronics response behavior of a Mu2e straw in
// several functions.
//
// $Id: StrawElectronics.cc,v 1.17 2014/09/22 12:23:28 brownd Exp $
// $Author: brownd $
// $Date: 2014/09/22 12:23:28 $
//
// Original author David Brown, LBNL
//
#include "TrackerConditions/inc/StrawElectronics.hh"
#include "GeneralUtilities/inc/DigitalFiltering.hh"
#include "cetlib_except/exception.h"
#include "TMath.h"
#include <math.h>
#include <algorithm>

#include <TFile.h>
#include <TH1F.h>

using namespace std;
namespace mu2e {
  using namespace TrkTypes;
  double StrawElectronics::_pC_per_uA_ns(1000.0); // unit conversion from pC/ns to microAmp.

  StrawElectronics::StrawElectronics(fhicl::ParameterSet const& pset) :
    _dVdI{pset.get<double>("thresholddVdI",1.8e4),
      pset.get<double>("adcdVdI",2.4e7),
      pset.get<double>("threshToAdc1dVdI",0),
      pset.get<double>("threshToAdc2dVdI",0) }, // mVolt/uAmps (transimpedance gain)
    _tdeadAnalog(pset.get<double>("DeadTimeAnalog",100.0)), // nsec dead after threshold crossing (pulse baseline restoration time)
    _tdeadDigital(pset.get<double>("DeadTimeDigital",100.0)), // nsec dead after threshold crossing (electronics processing time)
    _vsat(pset.get<double>("SaturationVoltage",120.0)), // mVolt
    _disp(pset.get<double>("Dispersion",1.0e-4)), // 0.1 ps/mm
    _vthresh(pset.get<double>("DiscriminatorThreshold",12.0)), //mVolt, post amplification
    _analognoise{pset.get<double>("thresholdAnalogNoise",3.0), //mVolt
      pset.get<double>("adcAnalogNoise",8.0),0,0},
    _ADCLSB(pset.get<double>("ADCLSB",0.3662)), //mVolt
    _maxADC(pset.get<int>("maxADC",4095)),
    _ADCped(pset.get<unsigned>("ADCPedestal",1393)),
    _nADC(pset.get<unsigned>("nADC",12)),
    _nADCpre(pset.get<unsigned>("nADCPresamples",4)),
    _ADCPeriod(pset.get<double>("ADCPeriod",20.0)), // nsec
    _ADCOffset(pset.get<double>("ADCOffset",2.0)), // nsec
    _TDCLSB(pset.get<double>("TDCLSB",0.015625)),  // nsec
    _maxTDC(pset.get<unsigned>("maxTDC",16777216)),
    _TOTLSB(pset.get<double>("TOTLSB",4.0)), //ns
    _maxTOT(pset.get<unsigned>("maxTOT",16)),
    _clockStart(pset.get<double>("clockStart",10.0)), // nsec
    _clockJitter(pset.get<double>("clockJitter",0.2)), // nsec
    _flashStart(pset.get<double>("FlashStart",0.0)), //nsec
    _flashEnd(pset.get<double>("FlashEnd",300.0)), // nsec
    _pmpEnergyScale(pset.get<double>("peakMinusPedestalEnergyScale",1.0)), // fudge factor for peak minus pedestal energy method

    _responseBins(pset.get<int>("ResponseBins",10000)),
    _sampleRate(pset.get<double>("SampleRate",1.0)), // ghz
    _saturationSampleFactor(pset.get<int>("SaturationSampleFactor",5)),
    _preampPoles(pset.get<vector<double> >("PreampPoles",vector<double>{160.,160.,6.0})),
    _preampZeros(pset.get<vector<double> >("PreampZeros",vector<double>{0.72343156})),
    _adcPoles(pset.get<vector<double> >("ADCPoles",vector<double>{6.24137023,4.6,30.0})),
    _adcZeros(pset.get<vector<double> >("ADCZeros",vector<double>{0.72343156})),
    _preampToAdc1Poles(pset.get<vector<double> >("PreampToAdc1Poles",vector<double>{6.24137032,30.0})),
    _preampToAdc1Zeros(pset.get<vector<double> >("PreampToAdc1Zeros",vector<double>{0.72343156})),
    _preampToAdc2Poles(pset.get<vector<double> >("PreampToAdc2Poles",vector<double>{4.6})),
    _preampToAdc2Zeros(pset.get<vector<double> >("PreampToAdc2Zeros",vector<double>{})),
    _ionNormalization(pset.get<double>("IonNormalization",1.0)),
    _ionSigma(pset.get<double>("IonSigma",2.0)),
    _ionT0(pset.get<double>("IonT0",6.0))
 {
    // precompute ion drift current pulses
    _currentPulse = std::vector<double>(_responseBins,0);
    _currentImpulse = std::vector<double>(_responseBins,0);
    _currentImpulse[_responseBins/2] = 1;
    double integral = 0;
    for (int i=0;i<_responseBins;i++){
      double t_gaus = (i-_responseBins/2)/_sampleRate;
      double val_gaus = 1/sqrt(TMath::TwoPi()*_ionSigma*_ionSigma)*exp(-(t_gaus*t_gaus)/(2*_ionSigma*_ionSigma));
      for (int j=0;j<_responseBins-i;j++){
        double t_tail = j/_sampleRate;
        double val = val_gaus / (t_tail + _ionT0);
        _currentPulse[i+j] += val;
        integral += val;
      }
    }
    for (int i=0;i<_responseBins;i++){
      // correct for sampleRate so that calculateResponse peak is independent of it
      // this combined with pC_per_uA_ns is the unit transform from pC to uA
      // do not renormalize since changed shape can change normalization???
      // normalization here is folded into dVdI
      _currentPulse[i] *= _sampleRate / _pC_per_uA_ns;
    }
    
    // calculate parameters for transfer function
    _preampResponse = std::vector<double>(_responseBins,0);
    _adcResponse = std::vector<double>(_responseBins,0);
    _preampToAdc1Response = std::vector<double>(_responseBins,0);
    _preampToAdc2Response = std::vector<double>(_responseBins,0);
    calculateResponse(_preampPoles,_preampZeros,_currentPulse,_preampResponse,_dVdI[thresh]);
    calculateResponse(_adcPoles,_adcZeros,_currentPulse,_adcResponse,_dVdI[adc]);
    calculateResponse(_preampToAdc1Poles,_preampToAdc1Zeros,_currentPulse,_preampToAdc1Response,_dVdI[satadc1]);
    calculateResponse(_preampToAdc2Poles,_preampToAdc2Zeros,_currentImpulse,_preampToAdc2Response,_dVdI[satadc2]);

    // now set other parameters
    _tmax[thresh] = 0;
    _linmax[thresh] = 0;
    _tmax[adc] = 0;
    _linmax[adc] = 0;
    _ttrunc[thresh] = (_responseBins/2)/_sampleRate;
    _ttrunc[adc] = (_responseBins/2)/_sampleRate;
    for (int i=0;i<_responseBins;i++){
      if (_preampResponse[i] > _linmax[thresh]){
        _linmax[thresh] = _preampResponse[i];
        _tmax[thresh] = (-_responseBins/2 + i)/_sampleRate;
      }
      if (_adcResponse[i] > _linmax[adc]){
        _linmax[adc] = _adcResponse[i];
        _tmax[adc] = (-_responseBins/2 + i)/_sampleRate;
      }
    }
  }

  StrawElectronics::~StrawElectronics() {}

  void StrawElectronics::calculateResponse(std::vector<double> &poles, std::vector<double> &zeros, std::vector<double> &input, std::vector<double> &response, double dVdI) {
    std::vector<double> za;
    std::vector<double> pa;
    for (size_t i=0;i<poles.size();i++){
      if (poles[i] != 0)
        pa.push_back(poles[i]*-1*TMath::TwoPi());
    }
    for (size_t i=0;i<zeros.size();i++){
      if (zeros[i] != 0)
        za.push_back(zeros[i]*-1*TMath::TwoPi());
    }
    std::vector<double> a(za.size()+1,0);
    std::vector<double> b(pa.size()+1,0);
    std::vector<double> aprime(std::max(za.size()+1,pa.size()+1),0);
    std::vector<double> bprime(std::max(za.size()+1,pa.size()+1),0);
    DigitalFiltering::zpk2tf(b,a,za,pa); 
    DigitalFiltering::bilinear(bprime,aprime,a,b,_sampleRate*1000.);

    // calculate impulse response
    for (size_t i=0;i<static_cast<size_t>(_responseBins);i++){
      response[i] = 0;
      for (size_t j=0;j<bprime.size();j++){
        if (i >= j){
          response[i] += input[i-j]*bprime[j];
        }
      }
      for (size_t j=1;j<aprime.size();j++){
        if (i >= j){
          response[i] += -1*response[i-j]*aprime[j];
        }
      }
    }

    for (int i=0;i<_responseBins;i++)
      response[i] *= dVdI;
  }

  double StrawElectronics::linearResponse(Path ipath, double time, double charge) const {
    int index = time*_sampleRate + _responseBins/2.;
    if ( index >= _responseBins)
      index = _responseBins-1;
    if (index < 0)
      index = 0;
    if (ipath == thresh)
      return charge * _preampResponse[index];
    else if (ipath == adc)
      return charge * _adcResponse[index];
    else if (ipath == satadc1)
      return charge * _preampToAdc1Response[index];
    else
      return charge * _preampToAdc2Response[index] * _saturationSampleFactor;
  }
  
  double StrawElectronics::saturatedResponse(double vlin) const {
    if (vlin < _vsat)
      return vlin;
    else
      return _vsat;
  }

  unsigned short StrawElectronics::adcResponse(double mvolts) const {
    return min(static_cast<unsigned short>(max(static_cast<int>(floor(mvolts/_ADCLSB)+_ADCped),0)),_maxADC);
  }

  unsigned long StrawElectronics::tdcResponse(double time) const {
    // Offset to when the TDC clock starts
    return min(static_cast<unsigned long>(max(static_cast<int>(floor((time-_clockStart)/_TDCLSB)),0)),_maxTDC);
  }

  
void StrawElectronics::digitizeWaveform(ADCVoltages const& wf, ADCWaveform& adc) const{
    if(wf.size() != _nADC)
      throw cet::exception("SIM") 
	<< "mu2e::StrawElectronics: wrong number of voltages to digitize" 
	<< endl;
    adc.clear();
    adc.reserve(_nADC);
    for(auto iwf=wf.begin();iwf!=wf.end();++iwf)
      adc.push_back(adcResponse(*iwf));
  }

  void StrawElectronics::digitizeTimes(TDCTimes const& times,TDCValues& tdc) const {
    for(size_t itime=0;itime<2;++itime)
      tdc[itime] = tdcResponse(times[itime]);
  }

  void StrawElectronics::adcTimes(double time, ADCTimes& adctimes) const {
// clock has a fixed phase; Assume we digitize with a fixed delay relative to the leading edge
    adctimes.clear();
    adctimes.reserve(_nADC);
// find the phase immediately proceeding this time.  Subtract presamples
    size_t phase = std::max((int)0,int(ceil(time/_ADCPeriod))-(int)_nADCpre);
    for(unsigned itime=0;itime<_nADC;++itime){
      adctimes.push_back((phase+itime)*_ADCPeriod+_ADCOffset);
    }
  }
  
  bool StrawElectronics::combineEnds(double t1, double t2) const {
    // currently two clock ticks are allowed for coincidence time
    int clockTicks1 = static_cast<int>(floor(t1-_clockStart)/_ADCPeriod);
    int clockTicks2 = static_cast<int>(floor(t1-_clockStart)/_ADCPeriod);
    return abs(clockTicks1-clockTicks2) < 2;
  }

  void StrawElectronics::tdcTimes(TDCValues const& tdc, TDCTimes& times) const {
    for(size_t itime=0;itime<2;++itime)
    // add back the time when the clock started
      times[itime] = tdc[itime]*_TDCLSB+_clockStart;
  }
  
  double StrawElectronics::adcVoltage(unsigned short adcval) const {
    return (adcval-_ADCped)*_ADCLSB;
  }

  double StrawElectronics::adcCurrent(unsigned short adcval) const {
  // this includes the effects from normalization of the pulse shape
    return adcVoltage(adcval)/_dVdI[adc];
  }

  double StrawElectronics::mypow(double val,unsigned n) {
    switch ( n ) {
      case 1:
	return val; break;
      case 2:
	return val*val; break;
      case 3:
	return val*val*val; break;
      case 4:
	return val*val*val*val; break;
      case 5:
	return val*val*val*val*val; break;
      default:
	return 0.0; break;
    }
  }

}
