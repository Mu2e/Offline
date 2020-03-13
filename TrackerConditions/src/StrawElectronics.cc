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
#include <cmath>
#include <complex>
#include <algorithm>

using namespace std;
namespace mu2e {
  using namespace TrkTypes;

  void StrawElectronics::calculateResponse(std::vector<double> const& poles, 
		std::vector<double> const& zeros, std::vector<double> const& input, 
		std::vector<double> &response) {

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

    // calculate gain at 160 MHz
    std::complex<double> w (0,160 * TMath::TwoPi());
    std::complex<double> numerator,denominator;
    for (size_t i=0;i<b.size();i++)
      denominator += std::complex<double>(b[i],0)*std::pow(w,b.size()-i-1);
    for (size_t i=0;i<a.size();i++)
      numerator += std::complex<double>(a[i],0)*std::pow(w,a.size()-i-1);

    double gain_160 = std::abs(numerator/denominator);

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
      response[i] *= 1 / gain_160;
  }

  double StrawElectronics::linearResponse(Straw const& straw, Path ipath, double time, double charge, double distance, bool forsaturation) const {
    int index = time*_sampleRate + _responseBins/2.;
    if ( index >= _responseBins)
      index = _responseBins-1;
    if (index < 0)
      index = 0;

    double straw_length = 2*straw.halfLength();
    double reflection_time = _reflectionTimeShift + (2*straw_length-2*distance)/_reflectionVelocity;
    int index_refl = (time - reflection_time)*_sampleRate + _responseBins/2.;
    if (index_refl >= _responseBins)
      index_refl = _responseBins-1;
    if (index_refl < 0)
      index_refl = 0;

    double reflection_scale = _reflectionFrac * exp(-(2*straw_length-2*distance)/_reflectionALength);

    int  distIndex = 0;
    for (size_t i=1;i<_wPoints.size()-1;i++){
      if (distance < _wPoints[i]._distance)
        break;
      distIndex = i;
    }
    double distFrac = 1 - (distance - _wPoints[distIndex]._distance)/(_wPoints[distIndex+1]._distance - _wPoints[distIndex]._distance);
    double p0, p1;
    if (ipath == thresh){
      if (forsaturation){
        p0 = _wPoints[distIndex]._preampToAdc1Response[index]      + _wPoints[distIndex]._preampToAdc1Response[index_refl]*reflection_scale;
        p1 = _wPoints[distIndex + 1]._preampToAdc1Response[index]  + _wPoints[distIndex + 1]._preampToAdc1Response[index_refl]*reflection_scale;
      }else{
        p0 = _wPoints[distIndex]._preampResponse[index]      + _wPoints[distIndex]._preampResponse[index_refl]*reflection_scale;
        p1 = _wPoints[distIndex + 1]._preampResponse[index]  + _wPoints[distIndex + 1]._preampResponse[index_refl]*reflection_scale;
      }
    }else{
      p0 = _wPoints[distIndex]._adcResponse[index]      + _wPoints[distIndex]._adcResponse[index_refl]*reflection_scale;
      p1 = _wPoints[distIndex + 1]._adcResponse[index]  + _wPoints[distIndex + 1]._adcResponse[index_refl]*reflection_scale;
    }
    return charge * ( p0 * distFrac + p1 * (1 - distFrac)) * _dVdI[ipath][straw.id().getStraw()];
  }

  double StrawElectronics::adcImpulseResponse(StrawId sid, double time, double charge) const {
    int index = time*_sampleRate + _responseBins/2.;
    if ( index >= _responseBins)
      index = _responseBins-1;
    if (index < 0)
      index = 0;
    return charge * _preampToAdc2Response[index] * _saturationSampleFactor * _dVdI[adc][sid.getStraw()]/_dVdI[thresh][sid.getStraw()];
  }
 
  double StrawElectronics::saturatedResponse(double vlin) const {
    if (vlin < _vsat)
      return vlin;
    else
      return _vsat;
  }

  double StrawElectronics::maxResponseTime(StrawId sid, Path ipath,double distance) const {
    int  distIndex = 0;
    for (size_t i=1;i<_wPoints.size()-1;i++){
      if (distance < _wPoints[i]._distance)
        break;
      distIndex = i;
    }
    double distFrac = 1 - (distance - _wPoints[distIndex]._distance)/(_wPoints[distIndex+1]._distance - _wPoints[distIndex]._distance);
    double p0 = _wPoints[distIndex]._tmax[ipath][sid.getStraw()];
    double p1 = _wPoints[distIndex + 1]._tmax[ipath][sid.getStraw()];
    
    return p0 * distFrac + p1 * (1 - distFrac);
  }

  double StrawElectronics::maxLinearResponse(StrawId sid, Path ipath,double distance,double charge) const {
    int  distIndex = 0;
    for (size_t i=1;i<_wPoints.size()-1;i++){
      if (distance < _wPoints[i]._distance)
        break;
      distIndex = i;
    }
    double distFrac = 1 - (distance - _wPoints[distIndex]._distance)/(_wPoints[distIndex+1]._distance - _wPoints[distIndex]._distance);
    double p0 = _wPoints[distIndex]._linmax[ipath][sid.getStraw()];
    double p1 = _wPoints[distIndex + 1]._linmax[ipath][sid.getStraw()];
 
    return charge * (p0 * distFrac + p1 * (1 - distFrac)) * _dVdI[ipath][sid.getStraw()];
  }

  ADCValue StrawElectronics::adcResponse(StrawId sid, double mvolts) const {
    return min(static_cast<ADCValue>(max(static_cast<int>(floor(mvolts/_ADCLSB)+_ADCped[sid.getStraw()]),0)),_maxADC);
  }

  TDCValue StrawElectronics::tdcResponse(double time) const {
    // here time should already be relative to event window marker arrival at ROC.
    // we add electronics absolute delay before TDC actually digitizeds
    double time_from_ewm = time+_electronicsTimeDelay;
    return min(static_cast<TDCValue>(max(static_cast<int>(floor((time_from_ewm)/_TDCLSB)),0)),_maxTDC);
  }

  
  void StrawElectronics::digitizeWaveform(StrawId sid, ADCVoltages const& wf, ADCWaveform& adc) const{
    if(wf.size() != adc.size()){
      throw cet::exception("SIM") 
	<< "mu2e::StrawElectronics: wrong number of voltages to digitize" 
	<< endl;
    }
    for(size_t iadc=0;iadc<adc.size();++iadc){
      adc.at(iadc) = adcResponse(sid, wf[iadc]);
    }
  }

  bool StrawElectronics::digitizeAllTimes(TDCTimes const& times,double mbtime, TDCValues& tdcs) const {
    for(size_t itime=0;itime<2;++itime)
      tdcs[itime] = tdcResponse(times[itime]);
    // test these times against time wrapping.  This calculation should be done at construction or init,
    // FIXME!
    bool notwrap(true);
    for(auto tdc : tdcs){
      notwrap &= tdc > tdcResponse(_flashStart - mbtime) && tdc < _flashStartTDC;
    }
    return notwrap;
  }

  bool StrawElectronics::digitizeTimes(TDCTimes const& times,TDCValues& tdcs) const {
    for(size_t itime=0;itime<2;++itime)
      tdcs[itime] = tdcResponse(times[itime]);
    // test bothe these times against the flash blanking. 
    bool notflash(true);
    for(auto tdc : tdcs){
      notflash &= tdc > _flashEndTDC && tdc < _flashStartTDC;
    }
    return notflash;
  }

  void StrawElectronics::adcTimes(double time, ADCTimes& adctimes) const {
// clock has a fixed phase; Assume we digitize with a fixed delay relative to the leading edge
    adctimes.clear();
    adctimes.reserve(TrkTypes::NADC);
// find the phase immediately proceeding this time.  Subtract presamples
    size_t phase = std::max((int)0,int(ceil(time/_ADCPeriod))-(int)_nADCpre);
    for(unsigned itime=0;itime<TrkTypes::NADC;++itime){
      adctimes.push_back((phase+itime)*_ADCPeriod+_ADCOffset);
    }
  }
  
  bool StrawElectronics::combineEnds(double t1, double t2) const {
    // currently two clock ticks are allowed for coincidence time
    int clockTicks1 = static_cast<int>(floor((t1)/_ADCPeriod));
    int clockTicks2 = static_cast<int>(floor((t2)/_ADCPeriod));
    return (unsigned)abs(clockTicks1-clockTicks2) < _maxtsep;
  }

  void StrawElectronics::uncalibrateTimes(TrkTypes::TDCTimes &times, const StrawId &id) const {
    times[StrawEnd::hv] -= _timeOffsetPanel[id.getPanel()] + _timeOffsetStrawHV[id.getStraw()];
    times[StrawEnd::cal] -= _timeOffsetPanel[id.getPanel()] + _timeOffsetStrawCal[id.getStraw()];
  }

  double StrawElectronics::adcVoltage(StrawId sid, ADCValue adcval) const {
    return (adcval-_ADCped[sid.getStraw()])*_ADCLSB;
  }

  double StrawElectronics::adcCurrent(StrawId sid, ADCValue adcval) const {
  // this includes the effects from normalization of the pulse shape
    return adcVoltage(sid,adcval)/_dVdI[adc][sid.getStraw()];
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

  void StrawElectronics::print(std::ostream& os) const {
    os << endl << "StrawElectronics parameters: "  << std::endl;
    for(size_t i=0; i<_dVdI.size(); i++) {
      auto const& aa = _dVdI[i];
      string ss = string("dVdI[") + to_string(i) + string("]");
      printVector(os,ss,aa);
    }

    os << "ttrunc["<<_ttrunc.size()<<"] = ";
    for(auto x: _ttrunc) os << x << " " ;
    os << endl;

    os << "tdeadAnalog = " << _tdeadAnalog << endl;
    os << "tdeadDigital = " << _tdeadDigital << endl;
    os << "vsat = " << _vsat << endl;
    printVector(os,"vthresh",_vthresh);
    os << "snoise = " << _snoise << endl;

    os << "analognoise["<<_analognoise.size()<<"] = ";
    for(auto x: _analognoise) os << x << " " ;
    os << endl;

    os << "ADCLSB = " << _ADCLSB << endl;
    os << "maxADC = " << _maxADC << endl;

    os << "ADCped ("<<_ADCped.size()<<") = ";
    for(auto x: _ADCped) os << x << " " ;
    os << endl;

    os << "nADCpre = " << _nADCpre << endl;
    os << "ADCPeriod = " << _ADCPeriod << endl;
    os << "ADCOffset = " << _ADCOffset << endl;
    os << "maxtsep = " << _maxtsep << endl;
    os << "TCoince = " << _TCoince << endl;
    os << "maxTDC = " << _maxTDC << endl;
    os << "maxTOT = " << _maxTOT << endl;
    os << "tdcResolution = " << _tdcResolution << endl;
    os << "electronicsTimeDelay = " << _electronicsTimeDelay << endl;
    os << "ewMarkerROCJitter = " << _ewMarkerROCJitter << endl;
    os << "flashStart = " << _flashStart << endl;
    os << "flashEnd = " << _flashEnd << endl;
    os << "flashClockSpeed = " << _flashClockSpeed << endl;
    os << "flashStartTDC = " << _flashStartTDC << endl;
    os << "flashEndTDC = " << _flashEndTDC << endl;
    os << "responseBins = " << _responseBins << endl;
    os << "sampleRate = " << _sampleRate << endl;
    os << "saturationSampleFactor = " << _saturationSampleFactor << endl;
    printVector(os,"preampPoles",_preampPoles);          
    printVector(os,"preampZeros",_preampZeros);
    printVector(os,"adcPoles",_adcPoles);
    printVector(os,"adcZeros",_adcZeros);
    printVector(os,"preampToAdc1Poles",_preampToAdc1Poles);
    printVector(os,"preampToAdc1Zeros",_preampToAdc1Zeros);
    printVector(os,"preampToAdc2Poles",_preampToAdc2Poles);
    printVector(os,"preampToAdc2Zeros",_preampToAdc2Zeros);
    printVector(os,"wireDistances",_wireDistances);
    printVector(os,"currentMeans",_currentMeans);
    printVector(os,"currentNormalizations",_currentNormalizations);
    printVector(os,"currentSigmas",_currentSigmas);
    printVector(os,"currentT0s",_currentT0s);
    printVector(os,"currentImpulse",_currentImpulse);
    printVector(os,"preampToAdc2Response",_preampToAdc2Response);  
    if( _wPoints.size()>0) {
      size_t i = _wPoints.size()/2;
      os << "wPoints midpoint entry ("<<i<<") : " << endl;
      os << "   distance = " << _wPoints[i]._distance << endl;
      os << "   mean = " << _wPoints[i]._mean << endl;
      os << "   normalization = " << _wPoints[i]._normalization << endl;
      os << "   sigma = " << _wPoints[i]._sigma << endl;
      os << "   t0 = " << _wPoints[i]._t0 << endl;
      printVector(os,"   currentPulse",_wPoints[i]._currentPulse);
      printVector(os,"   preampResponse",_wPoints[i]._preampResponse);
      printVector(os,"   adcResponse",_wPoints[i]._adcResponse);
      printVector(os,"   preampToAdc1Response",_wPoints[i]._preampToAdc1Response);
    }
    os << "clusterLookbackTime = " << _clusterLookbackTime << endl;
    printVector(os,"timeOffsetPanel",_timeOffsetPanel);
    printVector(os,"timeOffsetStrawHV",_timeOffsetStrawHV);
    printVector(os,"timeOffsetStrawCal",_timeOffsetStrawCal);

  }

  void StrawElectronics::printVector(std::ostream& os, std::string const& name, 
				    std::vector<double> const& a) const {
    size_t n = a.size();
    if(n<=4) {
      os << name << " ("<<n<<") = ";
      for(auto x : a) os << x << " ";
      os << endl;
    } else {
      os << name <<" ("<<n<<") = " 
	 << a[0] << " " << a[1] << " ... " 
	 << a[n-2] << " " << a[n-1] << endl;
    }
  }

}
