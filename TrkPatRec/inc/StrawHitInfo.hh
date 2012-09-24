//
// $Id: StrawHitInfo.hh,v 1.3 2012/09/24 18:39:55 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/09/24 18:39:55 $
//
// struct for hit diagnostics
#ifndef StrawHitInfo_hh
#define StrawHitInfo_hh
//#include "BaBar/BaBar.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "Rtypes.h"

namespace mu2e {
  struct StrawHitInfo {
    CLHEP::Hep3Vector _pos;
    Float_t _edep, _time, _corrtime;
    Int_t _device, _sector, _layer, _straw;
    Bool_t _vloose, _loose, _tight, _delta;
    CLHEP::Hep3Vector _mcpos;
    Int_t _mcpdg, _mcgen, _mcproc;
    Float_t _mcedep, _mctime, _mct0, _mcmom, _mctd;
    Int_t _hflag;
    Float_t _hgd; // MVA output of generalized distance
    Float_t _dphi, _drho, _dt; // MVA inputs
// root
    ClassDef(StrawHitInfo,1)
  };
}
#endif
