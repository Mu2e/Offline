//
// StrawDigi is the offline representation of the raw data readout from a straw
//
// Original author David Brown, LBNL
//
#include "RecoDataProducts/inc/StrawDigi.hh"

namespace mu2e {
  using namespace TrkTypes;
  StrawDigi::StrawDigi() : _strawid(0)
  {
  }
  StrawDigi::StrawDigi(StrawId sid, TDCValues tdc, TOTValues tot, ADCWaveform const& adc) : _strawid(sid),
//  _tdc(tdc),
  _adc(adc)
  {
    for(size_t itdc=0;itdc<2;++itdc){
      _tdc[itdc] = tdc[itdc];
      _tot[itdc] = tot[itdc];
    }
  }
  
  StrawDigi::StrawDigi(StrawDigi const& other) : _strawid(other._strawid),
     _flag(other._flag),_adc(other._adc)
  {
    for(size_t itdc=0;itdc<2;++itdc){
      _tdc[itdc] = other._tdc[itdc];
      _tot[itdc] = other._tot[itdc];
    }
  }

  StrawDigi& StrawDigi::operator=(StrawDigi const& other) {
    if(this != &other){
      _strawid = other._strawid;
      _flag = other._flag;
      _adc = other._adc;
      for(size_t itdc=0;itdc<2;++itdc){
	_tdc[itdc] = other._tdc[itdc];
	_tot[itdc] = other._tot[itdc];
      }
    }
    return *this;
  }
}

