#ifndef TrkDiag_BkgHitInfo_HH
#define TrkDiag_BkgHitInfo_HH
#include "TrkDiag/inc/StrawHitInfo.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "Rtypes.h"
namespace mu2e {
// extend StrawHitInfo struct with additional info on background hits
  struct BkgHitInfo : public StrawHitInfo {
    Float_t _gdist; // generalized distance 
    Float_t _rrho; // transverse distance to cluster center
    Float_t _rerr; // estimated error on transverse distance
    Int_t _index; // index into StrawHit collection
    CLHEP::Hep3Vector _rpos; // relative position of hit WRT cluster center
    Bool_t _active; // hit is active in cluster
    Bool_t _cbkg; // hit is classified backg
  };
}
#endif
