//
// StrawElectronics collects the electronics response behavior of a Mu2e straw in
// several functions.
//
// $Id: StrawElectronics.cc,v 1.11 2014/01/16 21:03:22 brownd Exp $
// $Author: brownd $
// $Date: 2014/01/16 21:03:22 $
//
// Original author David Brown, LBNL
//
#include "TrackerMC/inc/StrawElectronics.hh"
#include "cetlib/exception.h"
#include <math.h>
#include <algorithm>

using namespace std;
namespace mu2e {

  StrawElectronics::StrawElectronics(fhicl::ParameterSet const& pset) :
    _dVdI(pset.get<double>("dVdI",200.0)), // mVolt/uAmps
    _tshape(pset.get<double>("ShapingTime",15.0)),
    _tpow(pset.get<double>("ShapingPower",1.0)),
    _tdead(pset.get<double>("DeadTime",10.0)), // nsec dead after threshold crossing (electronics processing time)
    _vmax(pset.get<double>("MaximumVoltage",1000.0)), // 1000 mV max
    _vsat(pset.get<double>("SaturationVoltage",800.0)),
    _disp(pset.get<double>("Dispersion",1.0e-4)), // 100 ps/m
    _vthresh(pset.get<double>("DiscriminatorThreshold",5.0)), //mVolt
    _vthreshnoise(pset.get<double>("DiscriminatorThresholdNoise",0.5)), //mVolt
    _ADCLSB(pset.get<double>("ADCLSB",0.25)),
    _maxADC(pset.get<unsigned>("maxADC",4095)),
    _ADCped(pset.get<unsigned>("ADCPedetal",64)),
    _nADC(pset.get<unsigned>("nADC",10)),
    _nADCpre(pset.get<unsigned>("nADCPresamples",2)),
    _ADCPeriod(pset.get<double>("ADCPeriod",20.0)),
    _ADCOffset(pset.get<double>("ADCOffset",2.0)),
    _TDCLSB(pset.get<double>("TDCLSB",0.037)),
    _maxTDC(pset.get<unsigned>("maxTDC",65535)),
    _maxDTDC(pset.get<unsigned>("maxDeltaTDC",250))
 {
    // calcluate normalization
    _teff = _tshape/_tpow;
    _norm = 1.0/(_tshape*tgamma(_tpow));
  }

  StrawElectronics::~StrawElectronics() {}

  double StrawElectronics::hitletResponse(double time,StrawHitlet const& hitlet) const {
// in future, split response by hitlet type. FIXME!!!
    double retval(0.0);
    // There is no response before the hitlets own time
    if(time >  hitlet.time()){
// response is relative to the hitlets own time, normalized by the shaping time.
      double dt = (time - hitlet.time())/_teff;
      retval = 1000.0*hitlet.charge()*_dVdI*_norm*pow(dt,_tpow)*exp(-dt);
// should make some smearing for dispersion here, based on the wire distance traveled. FIXME!!
    }
    return retval;
  }

  double StrawElectronics::saturatedResponse(double vlin) const {
    double retval(0.0);
  // up to saturation voltage the response is linear
    if(vlin < _vsat)
      retval = vlin;
    else
      retval = _vmax - (_vmax-_vsat)*exp(-(vlin-_vsat)/(_vmax-_vsat));
    return retval;
  }

  unsigned short StrawElectronics::adcResponse(double mvolts) const {
    static unsigned short zero(0);
    return min(max(static_cast<unsigned short>(floor(mvolts/_ADCLSB)+_ADCped),zero),_maxADC);
  }

  unsigned long StrawElectronics::tdcResponse(double time) const {
    static unsigned long zero(0);
    return min(max(static_cast<unsigned long>(floor(time/_TDCLSB)),zero),_maxTDC);
  }

  void StrawElectronics::digitizeWaveform(vector<double> const& wf,StrawDigi::ADCWaveform& adc) const{
    if(wf.size() != _nADC)      
      throw cet::exception("SIM") 
	<< "mu2e::StrawElectronics: wrong number of voltages to digitize" 
	<< endl;
    adc.clear();
    adc.reserve(_nADC);
    for(auto iwf=wf.begin();iwf!=wf.end();++iwf)
      adc.push_back(adcResponse(*iwf));
  }

  void StrawElectronics::digitizeTimes(array<double,2> const& times,StrawDigi::TDCValues& tdc) const {
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
    return fabs(t1-t2)/_TDCLSB < _maxDTDC;
  }

  void StrawElectronics::tdcTimes(StrawDigi::TDCValues const& tdc, std::array<double,2>& times) const {
    for(size_t itime=0;itime<2;++itime)
      times[itime] = tdc[itime]*_TDCLSB;
  }
  
  double StrawElectronics::adcVoltage(unsigned short adcval) const {
    return (adcval-_ADCped)*_ADCLSB;
  }

  double StrawElectronics::adcCurrent(unsigned short adcval) const {
  // this includes the effects from normalization of the pulse shape
    return adcVoltage(adcval)/_dVdI;
  }
}
