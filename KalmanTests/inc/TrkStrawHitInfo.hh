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
 
 // simple structs for diagnostics
  struct TrkStrawHitInfo {
    TrkStrawHitInfo();
    Int_t _active, _usable;
    Int_t _device, _sector, _layer, _straw;
    Int_t _nplane, _npanel, _nlayer;
    Float_t _z, _phi, _rho;
    Float_t _resid, _residerr, _rdrift, _rdrifterr, _trklen;
    Float_t _doca, _exerr, _penerr, _t0, _t0err;
    Float_t _ht, _tddist, _tdderr, _hlen, _wdot;
    Float_t _edep, _dx;
    Int_t _ambig;
  };

  struct TrkStrawHitInfoMC {
    TrkStrawHitInfoMC();
    Int_t _mcpdg, _mcgen, _mcproc, _mcrel;
    Float_t _mct0, _mcht, _mcdist, _mclen;
    Float_t _mcedep;
    Float_t _mcr, _mcphi;
    Int_t _mcambig;
    Bool_t _xtalk;
  };
 
 // deprecated struct, will go away soon
 struct TrkStrawHitInfo_old{
    Int_t _active, _usable, _device, _sector, _layer, _straw;
    Float_t _z, _phi, _rho;
    Float_t _resid, _residerr, _rdrift, _rdrifterr, _trklen;
    Float_t _doca, _exerr, _penerr, _t0, _t0err;
    Float_t _ht, _tddist, _tdderr, _hlen;
    Float_t _edep, _dx;
    Int_t _ambig;
    Int_t _mcn, _mcnunique, _mcppdg, _mcpgen, _mcpproc;
    Int_t _mcpdg, _mcgen, _mcproc;
    Float_t _mct0, _mcht, _mcdist, _mclen;
    Float_t _mcedep;
    Int_t _mcambig;
    Bool_t _xtalk;
  };
}
#endif
