#ifndef TrkDiag_ComboHitInfo_HH
#define TrkDiag_ComboHitInfo_HH
#include "Rtypes.h"
#include "Offline/DataProducts/inc/GenVector.hh"
namespace mu2e {
  // info about each hit in the combo hit
  struct ComboHitInfo {

    XYZVectorF _pos; // position of this hit
    XYVectorF _udir; // U direction of this hit
    Float_t _du; // distance from this hit to the parent ComboHit in U
    Float_t _dv; // distance from this hit to the parent ComboHit in V
    Float_t _dw; // distance from this hit to the parent ComboHit in W
    Float_t _uvar; // estimated variance along the U direction
    Float_t _vvar; // estimated variance along the V direction
    Float_t _wvar = std::numeric_limits<Float_t>::max(); // estimated variance along the W direction
    Float_t _tvar; // estimated time variance
    Float_t _thit = 0.0; // hit time
    Int_t _straw; // straw
    Int_t _upanel; // unique panel
    Int_t _nch; // hit counts
    Int_t _nsh;
  };
  struct ComboHitInfoMC {
    Int_t _rel; // relation to the 1st hit
    XYZVectorF _mcpos; // position of this hit
  };

}
#endif
