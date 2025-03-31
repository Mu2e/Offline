//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
#include "Offline/CalPatRec/inc/DeltaFinder_types.hh"

namespace mu2e {
  using  CalPatRec::HitData_t;

  float DeltaSeed::fSigT2 = 8*8;    // resolution in tCorr, squared, ns^2

//-----------------------------------------------------------------------------
// 'Hd1' could be nullptr - this is the case when a single hit is picked up
// in a station based on a prediction made from another station
// the hit could be overlapping with a hit associated with a found seed
// such 'picked-up' seeds could be identified ... how ?
// fIndex and fStation are assumed to be already set
//-----------------------------------------------------------------------------
  void DeltaSeed::Init(HitData_t* Hd0, HitData_t* Hd1, float Xc, float Yc, float Zc) {
    fSnx2             = 0;
    fSnxy             = 0;
    fSny2             = 0;
    fSnxr             = 0;
    fSnyr             = 0;

    fSumEDep          = 0;
    fSumT             = 0;
    fSumT2            = 0;

    fGood             =  1;

    int face0         = Hd0->fZFace;
    fSFace[0]         = face0;
    int face1         = (Hd1) ? Hd1->fZFace : -1;
    fSFace[1]         = face1;

    fDeltaIndex       = -1;
    fProtonIndex      = -1;

    for (int face=0; face<kNFaces; face++) {
      fFaceProcessed[face] = 0;
      fHitData      [face] = NULL;
    }

    fHitData[face0]         = Hd0;
    fChi21                  = Hd0->fChi2Min;

    if (Hd1) {
      fHitData[face1]       = Hd1;
      fType                 =  10*face0+face1;
      fNHits                =  2;
      fNStrawHits           = Hd0->fHit->nStrawHits()+Hd1->fHit->nStrawHits();
      fFaceProcessed[face0] = 1;
      fFaceProcessed[face1] = 1;
      fChi22                = Hd1->fChi2Min;
    }
    else {
      fType                 =  10*face0;  // may want to revisit
      fNHits                =  1;
      fNStrawHits           = Hd0->fHit->nStrawHits();
      fFaceProcessed[face0] =  1;
      fChi22                =  0;
    }

    fChi2Par          = fChi21+fChi22;
    fChi2Perp         = 0;
    fMinHitTime       = Hd0->fCorrTime;
    fMaxHitTime       = fMinHitTime;

    if (Hd1) {
      float ct2 = Hd1->fCorrTime;

      if (ct2 >= fMinHitTime) fMaxHitTime = ct2;
      else                    fMinHitTime = ct2;
    }
//-----------------------------------------------------------------------------
// coordinate information
//-----------------------------------------------------------------------------
    CofM.SetXYZ(Xc,Yc,Zc);

    for (int i=0; i<2; i++) {
      int face = fSFace[i];
      if (face < 0)                                                 continue;
      const HitData_t* hd = fHitData[face];
      const ComboHit*  ch = hd->fHit;
      if (hd) {
        fSnx2    += hd->fNx2;
        fSnxy    += hd->fNxy;
        fSny2    += hd->fNy2;
        fSnxr    += hd->fNxr;
        fSnyr    += hd->fNyr;

        fSumEDep += ch->energyDep()*ch->nStrawHits();
        fSumT    += hd->fCorrTime;
        fSumT2   += hd->fCorrTime*hd->fCorrTime;
      }
    }
  }

//-----------------------------------------------------------------------------
// add hit
//-----------------------------------------------------------------------------
    void DeltaSeed::AddHit(HitData_t* Hd) {

      int face             = Hd->fZFace;
      fHitData[face]       = Hd;
      fFaceProcessed[face] = 1;

      if (Hd->fCorrTime < fMinHitTime) fMinHitTime = Hd->fCorrTime;
      if (Hd->fCorrTime > fMaxHitTime) fMaxHitTime = Hd->fCorrTime;

      const ComboHit* ch = Hd->fHit;

      fNHits          += 1;
      fNStrawHits     += ch->nStrawHits();
//-----------------------------------------------------------------------------
// in parallel, update coordinate sums
//-----------------------------------------------------------------------------
      fSnx2    += Hd->fNx2;
      fSnxy    += Hd->fNxy;
      fSny2    += Hd->fNy2;
      fSnxr    += Hd->fNxr;
      fSnyr    += Hd->fNyr;

      fSumEDep += ch->energyDep()*ch->nStrawHits();
      fSumT    += Hd->fCorrTime;
      fSumT2   += Hd->fCorrTime*Hd->fCorrTime;
    }
//-----------------------------------------------------------------------------
// replace first hit with another one founce in the same face..
//-----------------------------------------------------------------------------
    void DeltaSeed::ReplaceFirstHit(HitData_t* Hd) {

      assert(fSFace[0] == Hd->fZFace);

      fHitData[fSFace[0]] = Hd;

      fChi21              = Hd->fChi2Min;
      fNStrawHits         = Hd->fHit->nStrawHits();
      fMinHitTime         = Hd->fCorrTime;
      fMaxHitTime         = fMinHitTime;

      fSumEDep            = Hd->fHit->energyDep()*Hd->fHit->nStrawHits();
      fSumT               = Hd->fCorrTime*Hd->fHit->nStrawHits();
      fSumT2              = Hd->fCorrTime*Hd->fCorrTime;
    }
//-----------------------------------------------------------------------------
// calculate Com and chi2's
//-----------------------------------------------------------------------------
  void DeltaSeed::CalculateCogAndChi2(float RCore, float SigmaR2) {
//-----------------------------------------------------------------------------
// update seed time and X and Y coordinates, accurate knowledge of Z is not very relevant
// if the seed has only two hits from initial intersection, there is no coordinates
// to redefine, chi2perp = 0 and chi2w is the sum of the two ...
//-----------------------------------------------------------------------------
      if ((fNHits == 2) and (fSFace[0] >= 0) and (fSFace[1] >= 0)) {
        fChi2Par  = fChi21+fChi22;
        fChi2Perp = 0;
        return;
      }
//-----------------------------------------------------------------------------
// the seed has more than two hits
//-----------------------------------------------------------------------------
      assert(fNHits > 2);

      double d  = fSnx2*fSny2-fSnxy*fSnxy;

      double xc = (fSnyr*fSnx2-fSnxr*fSnxy)/d;
      double yc = (fSnyr*fSnxy-fSnxr*fSny2)/d;

      CofM.SetX(xc);
      CofM.SetY(yc);
//-----------------------------------------------------------------------------
// calculate seed chi2 - can this be optimized ?
//-----------------------------------------------------------------------------
      fChi2Par  = 0;
      fChi2Perp = 0;

      for (int face=0; face<kNFaces; face++) {
        const HitData_t* hd = fHitData[face];
        if (hd) {
          double dx = hd->fX-xc;
          double dy = hd->fY-yc;
//-----------------------------------------------------------------------------
// split into wire parallel and perpendicular components
//-----------------------------------------------------------------------------
          float dxy_dot_w = dx*hd->fWx+dy*hd->fWy;
          float dxy_dot_n = dx*hd->fWy-dy*hd->fWx;

          float chi2_par  = (dxy_dot_w*dxy_dot_w)/(SigmaR2+hd->fSigW2);
          float drr       = fmax(fabs(dxy_dot_n)-RCore,0);
          float chi2_perp = (drr*drr)/SigmaR2;
          fChi2Par       += chi2_par;
          fChi2Perp      += chi2_perp;
        }
      }
    }
//-----------------------------------------------------------------------------
// utility: calculate chi2's, call for N>= 2 hit seeds
// if there are only two hits, Chi2Perp = 0
//-----------------------------------------------------------------------------
  void DeltaSeed::Chi2(float Xc, float Yc, float RCore, float SigmaR2, float& Chi2Par, float& Chi2Perp) {

    Chi2Par  = 0;
    Chi2Perp = 0;

    if ((fNHits == 2) and (fSFace[0] >= 0) and (fSFace[1] >= 0)) {
      fChi2Par  = fChi21+fChi22;
      return;
    }

   for (int face=0; face<kNFaces; face++) {
      const HitData_t* hd =fHitData[face];
      if (hd == nullptr)                                              continue;
//-----------------------------------------------------------------------------
// split chi^2 into parallel and perpendicular to the wire components
//-----------------------------------------------------------------------------
      float dx        = hd->fX-Xc;
      float dy        = hd->fY-Yc;

      float dxy_dot_w = dx*hd->fWx+dy*hd->fWy;
      float dxy_dot_n = dx*hd->fWy-dy*hd->fWx;

      float  chi2_par = (dxy_dot_w*dxy_dot_w)/(SigmaR2+hd->fSigW2);
      float drr       = fmax(fabs(dxy_dot_n)-RCore,0);
      float chi2_perp = (drr*drr)/SigmaR2;

      Chi2Par        += chi2_par;
      Chi2Perp       += chi2_perp;
    }
  }
}
