//
// StrawDigi is the offline representation of the raw data readout from a straw
// $Id: StrawDigi.cc,v 1.1 2013/12/07 19:50:42 brownd Exp $
// $Author: brownd $
// $Date: 2013/12/07 19:50:42 $
//
// Original author David Brown, LBNL
//
#include "RecoDataProducts/inc/StrawDigi.hh"

namespace mu2e {
  using namespace TrkTypes;
  StrawDigi::StrawDigi() : _strawIndex(0)
  {
  }
  StrawDigi::StrawDigi(StrawIndex index, TDCValues tdc, TOTValues tot, ADCWaveform const& adc) : _strawIndex(index),
//  _tdc(tdc),
  _adc(adc)
  {
    for(size_t itdc=0;itdc<2;++itdc){
      _tdc[itdc] = tdc[itdc];
      _tot[itdc] = tot[itdc];
    }
  }
  
  StrawDigi::StrawDigi(StrawDigi const& other) : _strawIndex(other._strawIndex),
  //  _tdc(other._tdc),
     _adc(other._adc)
  {
    for(size_t itdc=0;itdc<2;++itdc){
      _tdc[itdc] = other._tdc[itdc];
      _tot[itdc] = other._tot[itdc];
    }
  }

  StrawDigi& StrawDigi::operator=(StrawDigi const& other) {
    if(this != &other){
      _strawIndex = other._strawIndex;
      for(size_t itdc=0;itdc<2;++itdc){
	_tdc[itdc] = other._tdc[itdc];
	_tot[itdc] = other._tot[itdc];
      }
      _adc = other._adc;
    }
    return *this;
  }
}

