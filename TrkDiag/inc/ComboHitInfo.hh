#ifndef TrkDiag_ComboHitInfo_HH
#define TrkDiag_ComboHitInfo_HH
#include "Rtypes.h"
#include "Offline/DataProducts/inc/GenVector.hh"
namespace mu2e {
  // info about each hit in the combo hit
  struct ComboHitInfo {

    XYZVectorF _pos; // position of this hit
    XYZVectorF _udir; // U direction of this hit
    Float_t _du; // distance from this hit to the parent ComboHit in U
    Float_t _dv; // distance from this hit to the parent ComboHit in V
    Float_t _dw; // distance from this hit to the parent ComboHit in W
    Float_t _ures; // estimated error along the U direction
    Float_t _vres; // estimated error along the V direction
    Float_t _wres; // estimated error along the W direction
    Float_t _tres; // estimated time error
    Float_t _thit; // hit time
    Int_t _strawid; // strawid info
    Int_t _panelid; // panelid info
    Int_t _nch; // hit counts
    Int_t _nsh;
  };
  struct ComboHitInfoMC {
    Int_t _rel; // relation to the 1st hit
    XYZVectorF _mcpos; // position of this hit
  };

}
#endif
