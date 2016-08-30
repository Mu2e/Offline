
#ifndef TrkDiag_HitInfoMC_HH
#define TrkDiag_HitInfoMC_HH
#include "Rtypes.h"

namespace mu2e {
  struct HitInfoMC {
    Int_t _pdg, _gen, _proc, _rel;
    void reset() { _pdg = _gen = _proc = _rel = -1; }
  };
}
#endif
