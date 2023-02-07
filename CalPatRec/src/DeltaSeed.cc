//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
#include "Offline/CalPatRec/inc/DeltaFinder_types.hh"

namespace mu2e {
  using  DeltaFinderTypes::HitData_t;
//-----------------------------------------------------------------------------
// 'Hd1' could be nullptr - this is the case when a single hit is picked up
// in a station based on a prediction made from another station
// the hit could be overlapping with a hit associated with a found seed
// such 'picked-up' seeds could be identified ... how ?
//-----------------------------------------------------------------------------
  DeltaSeed::DeltaSeed(int Index, int Station, int Face0, HitData_t* Hd0, int Face1, HitData_t* Hd1) {
    Init(Index,Station,Face0,Hd0,Face1,Hd1);
  }


//-----------------------------------------------------------------------------
// initialization body
//-----------------------------------------------------------------------------
  void DeltaSeed::Init(int Index, int Station, int Face0, HitData_t* Hd0, int Face1, HitData_t* Hd1) {
    fSnx2             = 0;
    fSnxy             = 0;
    fSny2             = 0;
    fSnxr             = 0;
    fSnyr             = 0;

    fSumEDep          = 0;
    fSumT             = 0;

    fIndex            = Index;
    fStation          = Station;
    fGood             =  1;
    fSFace[0]         = Face0;
    fSFace[1]         = Face1;

    fDeltaIndex       = -1;
    fChi2Delta        = -1.;

    for (int face=0; face<kNFaces; face++) {
      fFaceProcessed[face] = 0;
      fHitData      [face] = NULL;
    }

    fHitData[Face0]         = Hd0;
    fChi21                  = Hd0->fChi2Min;

    if (Face1 >= 0) {
      fHitData[Face1]       = Hd1;
      fType                 =  10*Face0+Face1;
      fNHits                =  2;
      fNStrawHits           = Hd0->fHit->nStrawHits()+Hd1->fHit->nStrawHits();
      fFaceProcessed[Face0] = 1;
      fFaceProcessed[Face1] = 1;
      fChi22                = Hd1->fChi2Min;
    }
    else {
      fType                 =  10*Face0;  // may want to revisit
      fNHits                =  1;
      fNStrawHits           = Hd0->fHit->nStrawHits();
      fFaceProcessed[Face0] =  1;
      fChi22                =  0;
    }

    fChi2Par          = fChi21+fChi22;
    fChi2Perp         = 0;
    fMinHitTime       = Hd0->fHit->correctedTime();
    fMaxHitTime       = fMinHitTime;

    if (Hd1) {
      float ct2 = Hd1->fHit->correctedTime();

      if (ct2 >= fMinHitTime) fMaxHitTime = ct2;
      else                    fMinHitTime = ct2;
    }
//-----------------------------------------------------------------------------
// coordinate information
//-----------------------------------------------------------------------------
    for (int i=0; i<2; i++) {
      int face = fSFace[i];
      if (face < 0)                                                 continue;
      const HitData_t* hd = fHitData[face];
      const ComboHit*  ch = hd->fHit;
      if (hd) {
        double x0 = ch->pos().x();
        double y0 = ch->pos().y();
        double nx = ch->wdir().x();
        double ny = ch->wdir().y();

        // this should be the hit wdir - check...
        double nr = x0*ny-y0*nx;

        fSnx2 += nx*nx;
        fSnxy += nx*ny;
        fSny2 += ny*ny;
        fSnxr += nx*nr;
        fSnyr += ny*nr;

        fSumEDep += ch->energyDep()*ch->nStrawHits();
        fSumT    += ch->correctedTime()*ch->nStrawHits();
      }
    }
  }

//-----------------------------------------------------------------------------
// add hit
//-----------------------------------------------------------------------------
    void DeltaSeed::AddHit(HitData_t* Hd, int Face) {

      fHitData[Face]       = Hd;
      fFaceProcessed[Face] = 1;

      if (Hd->fCorrTime < fMinHitTime) fMinHitTime = Hd->fCorrTime;
      if (Hd->fCorrTime > fMaxHitTime) fMaxHitTime = Hd->fCorrTime;

      const ComboHit* ch = Hd->fHit;

      fNHits          += 1;
      fNStrawHits     += ch->nStrawHits();
//-----------------------------------------------------------------------------
// in parallel, update coordinate sums
//-----------------------------------------------------------------------------
      double x0        = ch->pos().x();
      double y0        = ch->pos().y();
      double nx        = ch->wdir().x();
      double ny        = ch->wdir().y();

      double nr        = x0*ny-y0*nx;

      fSnx2 += nx*nx;
      fSnxy += nx*ny;
      fSny2 += ny*ny;
      fSnxr += nx*nr;
      fSnyr += ny*nr;

      fSumEDep += ch->energyDep()*ch->nStrawHits();
      fSumT    += ch->correctedTime()*ch->nStrawHits();
    }
//-----------------------------------------------------------------------------
// replace first hit with another one founce in the same face..
//-----------------------------------------------------------------------------
    void DeltaSeed::ReplaceFirstHit(HitData_t* Hd) {

      fHitData[fSFace[0]] = Hd;

      fChi21              = Hd->fChi2Min;
      fNStrawHits         = Hd->fHit->nStrawHits();
      fMinHitTime         = Hd->fHit->correctedTime();
      fMaxHitTime         = fMinHitTime;

      fSumEDep            = Hd->fHit->energyDep()*Hd->fHit->nStrawHits();
      fSumT               = Hd->fHit->correctedTime()*Hd->fHit->nStrawHits();
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
        const HitData_t* hd =fHitData[face];
        if (hd) {
          double dx = hd->fHit->pos().x()-xc;
          double dy = hd->fHit->pos().y()-yc;
//-----------------------------------------------------------------------------
// split into wire parallel and perpendicular components
//-----------------------------------------------------------------------------
          const XYZVectorF& wdir = hd->fHit->wdir();

          float dxy_dot_w = dx*wdir.x()+dy*wdir.y();
          float dxy_dot_n = dx*wdir.y()-dy*wdir.x();

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
      float dx = hd->fHit->pos().x()-Xc;
      float dy = hd->fHit->pos().y()-Yc;
//-----------------------------------------------------------------------------
// split into wire parallel and perpendicular components
//-----------------------------------------------------------------------------
      const XYZVectorF& wdir = hd->fHit->wdir();

      float dxy_dot_w = dx*wdir.x()+dy*wdir.y();
      float dxy_dot_n = dx*wdir.y()-dy*wdir.x();

      float  chi2_par = (dxy_dot_w*dxy_dot_w)/(SigmaR2+hd->fSigW2);
      float drr       = fmax(fabs(dxy_dot_n)-RCore,0);
      float chi2_perp = (drr*drr)/SigmaR2;

      Chi2Par        += chi2_par;
      Chi2Perp       += chi2_perp;
    }
  }
}
