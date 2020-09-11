//
//
// struct for hit diagnostics
#ifndef StrawHitInfo_hh
#define StrawHitInfo_hh
//#include "BTrk/BaBar/BaBar.hh"
#include "DataProducts/inc/XYZVec.hh"
#include "Rtypes.h"

namespace mu2e {
  struct StrawHitInfo {
    XYZVec _pos;
    Float_t _edep, _time, _rho;
    Float_t _wres, _tres, _chisq, _stdt, _dist;
    Int_t _plane, _panel, _layer, _straw;
    Bool_t _bkg, _isolated, _tsel, _esel, _rsel, _stereo, _tdiv, _strawxtalk, _elecxtalk;
    Int_t _relation;
    XYZVec _mcpos;
    Int_t _mcpdg, _mcgen, _mcproc, _mcid;
    Float_t _mcedep, _mctime, _mct0, _mcht, _mcmom, _mctd;
  };
}
#endif
