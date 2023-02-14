///////////////////////////////////////////////////////////////////////////////
// ComboHit.hh needs the definition of ProductID ...
///////////////////////////////////////////////////////////////////////////////
#ifndef __CalPatRec_HitData_t_hh__
#define __CalPatRec_HitData_t_hh__

#include "canvas/Persistency/Provenance/ProductID.h"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"

namespace mu2e {
  class DeltaSeed;

  namespace DeltaFinderTypes {
    struct HitData_t {
      const ComboHit*         fHit;
      DeltaSeed*              fSeed;           // nullptr if not associated...
      int                     fUsed;           // TBD
      int                     fZFace;          // z-ordered face (for printing)
      int                     fDeltaIndex;     // is it really needed? - yes!
      float                   fChi2Min;
      float                   fSigW2;          // cached resolution^2 along the wire
      float                   fCorrTime;       // cached hit corrected time
      float                   fX;
      float                   fY;
      float                   fWx;
      float                   fWy;
      float                   fNr;
      float                   fNx2;
      float                   fNxy;
      float                   fNy2;
      float                   fNxr;
      float                   fNyr;

      HitData_t(const ComboHit* Hit,int ZFace) {
        fHit         = Hit;
        fSeed        = nullptr;
        fUsed        = 0;
        fZFace       = ZFace;
        fChi2Min     = 99999.0;
        float sigw   =  Hit->posRes(ComboHit::wire);
        fSigW2       = sigw*sigw;
        fCorrTime    = Hit->correctedTime();
        fX           = Hit->pos ().x();
        fY           = Hit->pos ().y();
        fWx          = Hit->wdir().x();
        fWy          = Hit->wdir().y();
        fNr          = fX*fWy-fY*fWx;
        fNx2         = fWx*fWx;
        fNxy         = fWx*fWy;
        fNy2         = fWy*fWy;
        fNxr         = fWx*fNr;
        fNyr         = fWy*fNr;
        fDeltaIndex  = -1;
      }

      int Used      () const { return fUsed      ; }
      int DeltaIndex() const { return fDeltaIndex; }
      int panelID   () const {
        printf("HitData_t::panelID : write me!\n");
        return -1;
      }
    };
  }
}
#endif
