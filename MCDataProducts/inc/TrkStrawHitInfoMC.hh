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
}
#endif
