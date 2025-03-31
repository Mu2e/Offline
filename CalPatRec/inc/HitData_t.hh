///////////////////////////////////////////////////////////////////////////////
// ComboHit.hh needs the definition of ProductID ...
///////////////////////////////////////////////////////////////////////////////
#ifndef CalPatRec_HitData_t_hh
#define CalPatRec_HitData_t_hh

#include <cmath>

#include "Offline/RecoDataProducts/inc/ComboHit.hh"

namespace mu2e {
  class DeltaSeed;

namespace CalPatRec {
  struct HitData_t {
    const mu2e::ComboHit*   fHit = nullptr;
    mu2e::DeltaSeed*        fSeed = nullptr; // nullptr if not associated...
    int                     fUsed = 0;       // TBD
    int                     fZFace = 0;      // z-ordered face (for printing)
    int                     fDeltaIndex = 0; // is it really needed? - yes!
    int                     fProtonIndex = 0;//
    float                   fChi2Min = 0.f;
    float                   fSigW2 = 0.f;    // cached resolution^2 along the wire
    float                   fCorrTime = 0.f; // cached hit corrected time
    float                   fX = 0.f;
    float                   fY = 0.f;
    float                   fWx = 0.f;
    float                   fWy = 0.f;
    float                   fNr = 0.f;
    float                   fNx2 = 0.f;
    float                   fNxy = 0.f;
    float                   fNy2 = 0.f;
    float                   fNxr = 0.f;
    float                   fNyr = 0.f;


    HitData_t(const mu2e::ComboHit* Hit,int ZFace) {
        fHit         = Hit;
        fSeed        = nullptr;
        fUsed        = 0;
        fZFace       = ZFace;
        fChi2Min     = 99999.0;
        float sigw   =  Hit->posRes(mu2e::ComboHit::wire);
        fSigW2       = sigw*sigw;
        fCorrTime    = Hit->correctedTime();
        fX           = Hit->pos ().x();
        fY           = Hit->pos ().y();
        fWx          = Hit->uDir2D().x();
        fWy          = Hit->uDir2D().y();
        fNr          = fX*fWy-fY*fWx;
        fNx2         = fWx*fWx;
        fNxy         = fWx*fWy;
        fNy2         = fWy*fWy;
        fNxr         = fWx*fNr;
        fNyr         = fWy*fNr;
        fDeltaIndex  = -1;
        fProtonIndex = -1;
      }

    int   Used       () const { return fUsed       ; }
    int   DeltaIndex () const { return fDeltaIndex ; }
    int   ProtonIndex() const { return fProtonIndex; }
    float Phi        () const { return atan2(fY,fX)  ; }

    void  setDeltaIndex (int Index) { fDeltaIndex  = Index; }
    void  setProtonIndex(int Index) { fProtonIndex = Index; }
    //-----------------------------------------------------------------------------
    // panel ID :
    //-----------------------------------------------------------------------------
    int   panelID    () const {
      int sid = fHit->strawId().asUint16();
      int panel_id = (sid >> 7) & 0x1ff;
      return panel_id;
    }
  };
}
}
#endif
