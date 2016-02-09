//
// Struct to hold pointers to MC data in event.
// $Id: TrkStrawHitInfo.hh,v 1.2 2014/09/22 12:13:17 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/09/22 12:13:17 $
//
#ifndef TrkStrawHitInfo_HH
#define TrkStrawHitInfo_HH
#include "Rtypes.h"
// simple structs
namespace mu2e 
{
 // simple struct for diagnostics
  struct TrkStrawHitInfo {
    Bool_t _active, _dhit, _dactive;
    Int_t _plane, _panel, _layer, _straw;
    Int_t _ambig;
    Float_t _z, _phi, _rho;
    Float_t _resid, _residerr, _rdrift, _rdrifterr, _trklen;
    Float_t _doca, _exerr, _penerr, _t0, _t0err;
    Float_t _ht, _tddist, _tdderr, _hlen, _wdot;
    Float_t _edep, _dx;
    TrkStrawHitInfo() : _active(false), _plane(-1),
    _panel(-1), _layer(-1), _straw(-1), _ambig(-1),
    _z(-1000.0), _phi(-1000.0), _rho(-1000.0),
    _resid(-1000.0), _residerr(-1000.0), _rdrift(-1000.0), _rdrifterr(-1000.0),
    _trklen(-1000.0),_doca(-1000.0), _exerr(-1000.0), _penerr(-1000.0),
    _t0(-1000.0), _t0err(-1000.0), _ht(-1000.0), _tddist(-1000.0), _tdderr(-1000.0),
    _hlen(-1000.0), _edep(-1000.0), _dx(-1000.0)  {}
  };
}
#endif
