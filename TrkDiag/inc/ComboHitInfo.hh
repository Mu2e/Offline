#ifndef TrkDiag_ComboHitInfo_HH
#define TrkDiag_ComboHitInfo_HH
#include "Rtypes.h"
#include "DataProducts/inc/XYZVec.hh"
namespace mu2e {
  // info about each hit in the combo hit
  struct ComboHitInfo {
    
    XYZVec _pos; // position of this hit
    XYZVec _wdir; // direction of this hit
    Float_t _wdist; // distance from wire center along this direction
    Float_t _wres; // estimated error along this direction
    Float_t _tres; // estimated error perpendicular to this direction
    Float_t _tdrift; // drift time estimate
    Float_t _thit; // hit time
    Int_t _strawid; // strawid info
    Int_t _panelid; // panelid info
    Float_t _dperp, _dwire; // distance from average to this hit
    Float_t _dterr; // estimated error on distance to center
    Float_t _dwerr; // estimated error along wire direction
    Float_t _dtime; // time diference to average
    Float_t _dedep; // energy dep difference
    Int_t _ds; // straw difference
    Int_t _dp; // panel difference
    Int_t _nch;
    Int_t _nsh;
  };
  struct ComboHitInfoMC {
    Int_t _rel; // relation to the 1st hit
    XYZVec _mcpos; // position of this hit
  };

}
#endif
