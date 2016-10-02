#ifndef TrkDiag_HitInfoMC_HH
#define TrkDiag_HitInfoMC_HH
#include "Rtypes.h"

namespace mu2e {
  struct HitInfoMC {
    Int_t _pdg, _gen, _proc, _rel;
    Float_t _t0;
    void reset() { _pdg = _gen = _proc = _rel = -1; _t0= 0.0;}
  };
}
#endif
