//
// $Id: StrawHitInfo.hh,v 1.4 2013/03/08 04:33:26 brownd Exp $
// $Author: brownd $ 
// $Date: 2013/03/08 04:33:26 $
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
    Float_t _edep, _time, _corrtime,_rho;
    Float_t _pres, _rres, _chisq;
    Int_t _device, _sector, _layer, _straw;
    Bool_t _esel, _rsel, _delta, _stereo;
    Int_t _relation;
    Bool_t _primary;
    CLHEP::Hep3Vector _mcpos;
    Int_t _mcpdg, _mcgen, _mcproc, _mcid;
    Float_t _mcedep, _mctime, _mct0, _mcmom, _mctd;
    Int_t _hflag;
    Float_t _hgd; // MVA output of generalized distance
    Float_t _dphi, _drho, _dt; // MVA inputs
// root
    ClassDef(StrawHitInfo,1)
  };
}
#endif
