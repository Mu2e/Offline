//
// StrawElectronics collects the electronics response behavior of a Mu2e straw in
// several functions.
//
// $Id: StrawElectronics.cc,v 1.2 2013/12/08 21:10:12 brownd Exp $
// $Author: brownd $
// $Date: 2013/12/08 21:10:12 $
//
// Original author David Brown, LBNL
//
#include "TrackerMC/inc/StrawElectronics.hh"
#include "cetlib/exception.h"

using namespace std;
namespace mu2e {

  StrawElectronics::StrawElectronics(fhicl::ParameterSet const& pset) :
    _dVdQ(pset.get<double>("dVdQ",1.0)),
    _trise(pset.get<double>("RiseTime",0.0)),
    _tfall(pset.get<double>("FallTime",50.0)),
    _vmax(pset.get<double>("MaximumVoltage",1000.0)), // 1000 mV max
    _vsat(pset.get<double>("SaturationVoltage",800.0)),
    _ADCLSB(pset.get<double>("ADCLSB",0.1)),
    _maxADC(pset.get<unsigned>("maxADC",1047)),
    _nADC(pset.get<unsigned>("nADC",10)),
    _ADCPeriod(pset.get<double>("ADCPeriod",15.4)),
    _ADCOffset(pset.get<double>("ADCOffset",40.0)),
    _TDCLSB(pset.get<double>("TDCLSB",0.037)),
    _maxTDC(pset.get<unsigned>("maxTDC",65535)),
    _maxDTDC(pset.get<unsigned>("maxDeltaTDC",200))
  {
  // insure times are positive
      if(_trise < 0.0 || _tfall < 0.0){
      throw cet::exception("SIM") 
	<< "mu2e::StrawElectronics: negative rise or fall time!" 
	<< endl;
      }
    // normalization is given by rise and fall times
    _norm = (_trise+_tfall)/(_tfall*_tfall);
    // relative time at maximum
    if(_trise > 0.0){
      _tmax = _trise*(log(_tfall + _trise) - log(_trise));
    } else
      _tmax = 0.0;
  }

  double StrawElectronics::hitletResponse(double time,StrawHitlet const& hitlet) const {
// in future, split response by hitlet type. FIXME!!!
    double retval(0.0);
    // There is no response before the hitlets own time
    if(time >  hitlet.time()){
// response is relative to the hitlets own time.
      double dt = time - hitlet.time();
   // 0 risetime is allowed
      double rise(1.0);
      if(_trise > 0.0)
	rise = (1.0-exp(-dt/_trise));
      double fall = exp(-dt/_tfall);
      retval = hitlet.charge()*_dVdQ*rise*fall;
// smear this by electronics noise FIXME!!!
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
    return min(max(static_cast<unsigned short>(floor(mvolts*_ADCLSB)),zero),_maxADC);
  }

  unsigned long StrawElectronics::tdcResponse(double time) const {
    static unsigned long zero(0);
    return min(max(static_cast<unsigned long>(floor(time*_TDCLSB)),zero),_maxTDC);
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
    adctimes.clear();
    adctimes.reserve(_nADC);
    for(unsigned itime=0;itime<_nADC;++itime){
      adctimes.push_back(time+_ADCOffset+itime*_ADCPeriod);
    }
  }
  
  bool StrawElectronics::combineEnds(double t1, double t2) const {
    return fabs(t1-t2)*_TDCLSB < _maxDTDC;
  }



}
