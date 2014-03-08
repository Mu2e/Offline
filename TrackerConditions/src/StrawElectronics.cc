//
// StrawElectronics collects the electronics response behavior of a Mu2e straw in
// several functions.
//
// $Id: StrawElectronics.cc,v 1.11 2014/03/08 00:55:21 brownd Exp $
// $Author: brownd $
// $Date: 2014/03/08 00:55:21 $
//
// Original author David Brown, LBNL
//
#include "TrackerConditions/inc/StrawElectronics.hh"
#include "cetlib/exception.h"
#include <math.h>
#include <algorithm>

using namespace std;
namespace mu2e {

  StrawElectronics::StrawElectronics(fhicl::ParameterSet const& pset) :
    _dVdI{pset.get<double>("thresholddVdI",85.0),
      pset.get<double>("adcdVdI",600.0) }, // mVolt/uAmps (transimpedance gain)
    _tshape{pset.get<double>("thresholdShapingTime",4.0),
      pset.get<double>("adcShapingTime",25.0) }, // nsec
    _tdead(pset.get<double>("DeadTime",60.0)), // nsec dead after threshold crossing (electronics processing time)
    _vmax(pset.get<double>("MaximumVoltage",1000.0)), // 1000 mVolt
    _vsat(pset.get<double>("SaturationVoltage",800.0)), // mVolt
    _disp(pset.get<double>("Dispersion",1.0e-4)), // 0.1 ps/mm
    _vthresh(pset.get<double>("DiscriminatorThreshold",20.0)), //mVolt, post amplification
    _vthreshnoise(pset.get<double>("DiscriminatorThresholdNoise",3.0)), //mVolt
    _ADCLSB(pset.get<double>("ADCLSB",0.25)), //mVolt
    _maxADC(pset.get<int>("maxADC",4095)),
    _ADCped(pset.get<unsigned>("ADCPedestal",128)),
    _nADC(pset.get<unsigned>("nADC",7)),
    _nADCpre(pset.get<unsigned>("nADCPresamples",2)),
    _ADCPeriod(pset.get<double>("ADCPeriod",20.0)), // nsec
    _ADCOffset(pset.get<double>("ADCOffset",2.0)), // nsec
    _TDCLSB(pset.get<double>("TDCLSB",0.037)),  // nsec
    _maxTDC(pset.get<unsigned>("maxTDC",65535)),
    _clockStart(pset.get<double>("clockStart",200.0)), // nsec
    _clockJitter(pset.get<double>("clockJitter",0.2)), // nsec
    _flashStart(pset.get<double>("FlashStart",0.0)), //nsec
    _flashEnd(pset.get<double>("FlashEnd",300.0)) // nsec
 {
    // calcluate normalization
   for(int ipath=0;ipath<2;++ipath){ 
     _fshape[ipath] = 1.0/_tshape[ipath];
     _norm[ipath] = 1000.0*_dVdI[ipath]/(_tshape[ipath]); // this includes unit conversion to microamps from charge in picoC and time in nsec
     _tmax[ipath] = _tshape[ipath]*pset.get<double>("MaxTime",10.0);
     _linmax[ipath] = linearResponse(static_cast<path>(ipath),_tshape[ipath],1.0); // response to unit charge (without saturation!)
   }
   _vdiff = _vmax-_vsat;
 }

  StrawElectronics::~StrawElectronics() {}

  double StrawElectronics::linearResponse(path ipath,double time,double charge) const {
    double retval(0.0);
    // There is no response before the hitlets own time
    if(time >  0.0 && time < _tmax[ipath]){
// response is relative to the time normalized by the shaping time.
      double tau = time*_fshape[ipath];
      retval = charge*_norm[ipath]*tau*exp(-tau);     
    }
    return retval;
  }

  double StrawElectronics::saturatedResponse(double vlin) const {
    double retval(0.0);
  // up to saturation voltage the response is linear
    if(vlin < _vsat)
      retval = vlin;
    else if (vlin > _vmax)
      retval = _vmax;
    else {
      retval = _vmax - _vdiff*exp(-(vlin-_vsat)/_vdiff);
    }
    return retval;
  }

  unsigned short StrawElectronics::adcResponse(double mvolts) const {
    return min(static_cast<unsigned short>(max(static_cast<int>(floor(mvolts/_ADCLSB)+_ADCped),0)),_maxADC);
  }

  unsigned long StrawElectronics::tdcResponse(double time) const {
    // Offset to when the TDC clock starts
    return min(static_cast<unsigned long>(max(static_cast<int>(floor((time-_clockStart)/_TDCLSB)),0)),_maxTDC);
  }
void StrawElectronics::digitizeWaveform(vector<double> const& wf, ADCWaveform& adc) const{
    if(wf.size() != _nADC)      
      throw cet::exception("SIM") 
	<< "mu2e::StrawElectronics: wrong number of voltages to digitize" 
	<< endl;
    adc.clear();
    adc.reserve(_nADC);
    for(auto iwf=wf.begin();iwf!=wf.end();++iwf)
      adc.push_back(adcResponse(*iwf));
  }

  void StrawElectronics::digitizeTimes(array<double,2> const& times,TDCValues& tdc) const {
    for(size_t itime=0;itime<2;++itime)
      tdc[itime] = tdcResponse(times[itime]);
  }

  void StrawElectronics::adcTimes(double time, vector<double>& adctimes) const {
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
  // use deadtime to define when ends should be combined: that
  // must be longer than the propagaion + resolution time!
    return fabs(t1-t2) < _tdead;
  }

  void StrawElectronics::tdcTimes(TDCValues const& tdc, std::array<double,2>& times) const {
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
