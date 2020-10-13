//
// StrawDigi is the offline representation of the raw data readout from a straw
//
// Original author David Brown, LBNL
//
#include "RecoDataProducts/inc/StrawDigi.hh"

namespace mu2e {
  StrawDigi::StrawDigi(StrawId sid, TrkTypes::TDCValues tdc, TrkTypes::TOTValues tot, TrkTypes::ADCValue pmp) : 
        _strawid(sid), _tdc(tdc), _tot(tot), _pmp(pmp) {};
};
