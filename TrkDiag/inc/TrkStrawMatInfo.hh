//
// Struct to describe straw materials on the Kalman fit
// $Id: TrkStrawMatInfo.hh,v 1.2 2014/09/22 12:13:17 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/09/22 12:13:17 $
//
#ifndef TrkStrawMatInfo_HH
#define TrkStrawMatInfo_HH
#include "Rtypes.h"
// simple structs
namespace mu2e 
{
 // simple struct for diagnostics
  struct TrkStrawMatInfo {
    Bool_t _active, _thit, _thita;
    Int_t _plane, _panel, _layer, _straw;
    Float_t _doca, _tlen;
    Float_t _dp, _radlen, _sigMS;
    TrkStrawMatInfo() : _active(false), _thit(false), _thita(false),
    _plane(-1), _panel(-1), _layer(-1), _straw(-1),
    _doca(-1000.0), _tlen(-1000.0),
    _dp(-1000.0), _radlen(-1000.0), _sigMS(-1000.0) {}
  };
}
#endif
