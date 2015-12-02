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
    Bool_t _active, _dhit, _dactive;
    Int_t _device, _panel, _layer, _straw;
    Int_t _nplane, _npanel, _nlayer;
    Int_t _ambig;
    Float_t _z, _phi, _rho;
    Float_t _resid, _residerr, _rdrift, _rdrifterr, _trklen;
    Float_t _doca, _exerr, _penerr, _t0, _t0err;
    Float_t _ht, _tddist, _tdderr, _hlen, _wdot;
    Float_t _edep, _dx;
    TrkStrawHitInfo() : _active(false), _device(-1),
    _panel(-1), _layer(-1), _straw(-1), _nplane(0), _npanel(0), _nlayer(0),_ambig(-1),
    _z(-1000.0), _phi(-1000.0), _rho(-1000.0),
    _resid(-1000.0), _residerr(-1000.0), _rdrift(-1000.0), _rdrifterr(-1000.0),
    _trklen(-1000.0),_doca(-1000.0), _exerr(-1000.0), _penerr(-1000.0),
    _t0(-1000.0), _t0err(-1000.0), _ht(-1000.0), _tddist(-1000.0), _tdderr(-1000.0),
    _hlen(-1000.0), _edep(-1000.0), _dx(-1000.0)  {}
  };

  struct TrkStrawHitInfoMC {
    Int_t _pdg, _gen, _proc, _rel;
    Float_t _t0, _ht, _dist, _len;
    Float_t _edep;
    Float_t _mom, _r, _phi;
    Int_t _ambig;
    Bool_t _xtalk;
    TrkStrawHitInfoMC() : _pdg(-1), _gen(-1), _proc(-1), 
    _t0(-1000.0), _ht(-1000.0), _dist(-1000.0), _len(-1000.0),
    _edep(-1000.0),_r(-1000.0),_phi(-1000.0),
    _ambig(-100), _xtalk(false) {}
  };
 
 // deprecated struct, will go away soon
 struct TrkStrawHitInfo_old{
    Int_t _active, _usable, _device, _panel, _layer, _straw;
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
