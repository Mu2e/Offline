#ifndef TrkDiag_HelixHitInfo_HH
#define TrkDiag_HelixHitInfo_HH
#include "TrkDiag/inc/HitInfoMC.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "Rtypes.h"
namespace mu2e {
  struct HelixHitInfo {
    Bool_t _outlier, _stereo, _tdiv, _resphi, _delta, _esel;
    CLHEP::Hep3Vector _hhpos, _hpos;
    Float_t _hhphi, _hphi;
    Float_t _werr, _terr, _dt;
    Float_t _whdot, _hrho;
    Float_t _dwire, _dtrans, _wres, _wtres, _chisq;
  };

  struct HelixHitInfoMC : public HitInfoMC {
    Float_t _hphi;
    Float_t _dwire, _dtrans;
    CLHEP::Hep3Vector _hpos;
  };
}
#endif
