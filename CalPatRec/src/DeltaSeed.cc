//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
#include "Offline/CalPatRec/inc/DeltaFinder_types.hh"
ClassImp(mu2e::DeltaSeed)

namespace mu2e {
  //  namespace DeltaFinderTypes {

//-----------------------------------------------------------------------------
    // DeltaSeed::DeltaSeed() : TObject() {
    //   printf("ERROR: DeltaSeed::DeltaSeed() should not be called\n");
    // }

//-----------------------------------------------------------------------------
    DeltaSeed::DeltaSeed(int Index, int Station, int Face0, HitData_t* Hd0, int Face1, HitData_t* Hd1):
      TObject(),
      fSnx2(0), fSnxy(0), fSny2(0),fSnxr(0),fSnyr(0),
      fSumEDep(0), fSumT(0)
    {
      fIndex            = Index;
      fStation          = Station;
      fGood             =  1;
      fSFace[0]         = Face0;
      fSFace[1]         = Face1;
      fChi2Delta        = -1.;

      for (int face=0; face<kNFaces; face++) {
        fFaceProcessed[face] = 0;
        hitData       [face] = NULL;
        fMcPart       [face] = NULL;
      }

      hitData[Face0]    = Hd0;
      fChi21            = Hd0->fChi2Min;

      if (Face1 >= 0) {
        hitData[Face1]        = Hd1;
        fType                 =  10*Face0+Face1;
        fNHits                =  2;
        fNStrawHits           = Hd0->fHit->nStrawHits()+Hd1->fHit->nStrawHits();
        fNFacesWithHits       =  2;
        fFaceProcessed[Face0] = 1;
        fFaceProcessed[Face1] = 1;
        fChi22                = Hd1->fChi2Min;
      }
      else {
        fType                 =  10*Face0;  // may want to revisit
        fNHits                =  1;
        fNStrawHits           = Hd0->fHit->nStrawHits();
        fNFacesWithHits       =  1;
        fFaceProcessed[Face0] =  1;
        fChi22                = -1;
      }

      fNHitsCE          =  0;
                                        // these ones used by teh diag tool
      fPreSeedMcPart[0] = nullptr;
      fPreSeedMcPart[1] = nullptr;

      fMinHitTime = Hd0->fHit->correctedTime();
      fMaxHitTime = fMinHitTime;

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
        const HitData_t* hd = hitData[face];
        const ComboHit* ch = hd->fHit;
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
//-----------------------------------------------------------------------------
// do't want a chi2 selection to do anything if the chi2 has not been calculated
//-----------------------------------------------------------------------------
      fDeltaIndex       = -1;
      fChi2All          = -1;
      fChi2Perp         = -1;
      fNSim             = -1;
      fNMom             = -1;
    }

//-----------------------------------------------------------------------------
// add hit
//-----------------------------------------------------------------------------
    void DeltaSeed::AddHit(HitData_t* Hd, int Face) {

      hitData[Face]        = Hd;
      fFaceProcessed[Face] = 1;

      if (Hd->fCorrTime < fMinHitTime) fMinHitTime = Hd->fCorrTime;
      if (Hd->fCorrTime > fMaxHitTime) fMaxHitTime = Hd->fCorrTime;

      const ComboHit* ch = Hd->fHit;

      fNFacesWithHits += 1;
      fNHits          += 1;
      fNStrawHits     += ch->nStrawHits();
//-----------------------------------------------------------------------------
// in parallel, update coordinate sums
//-----------------------------------------------------------------------------
      double x0                = ch->pos().x();
      double y0                = ch->pos().y();
      double nx                = ch->wdir().x();
      double ny                = ch->wdir().y();

      double nr                = x0*ny-y0*nx;

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

      hitData[fSFace[0]] = Hd;

      fChi21             = Hd->fChi2Min;
      fNStrawHits        = Hd->fHit->nStrawHits();
      fMinHitTime        = Hd->fHit->correctedTime();
      fMaxHitTime        = fMinHitTime;

      fSumEDep           = Hd->fHit->energyDep()*Hd->fHit->nStrawHits();
      fSumT              = Hd->fHit->correctedTime()*Hd->fHit->nStrawHits();
    }
//-----------------------------------------------------------------------------
// calculate Com and chi2's
//-----------------------------------------------------------------------------
    void DeltaSeed::CalculateCogAndChi2(double SigmaR2) {
//-----------------------------------------------------------------------------
// update seed time and X and Y coordinates, accurate knowledge of Z is not very relevant
// if the seed has only two hits from initial intersection, there is no coordinates
// to redefine, chi2perp = 0 and chi2w is the sum of the two ...
//-----------------------------------------------------------------------------
      if ((fNHits == 2) and (fSFace[0] >= 0) and (fSFace[1] >= 0)) {
        fChi2All  = fChi21+fChi22;
        fChi2Perp = 0;
        return;
      }
//-----------------------------------------------------------------------------
// the seed has more than two hits
//-----------------------------------------------------------------------------
      assert(fNHits > 2);

      double nxym, nx2m, ny2m, nxrm, nyrm;

      nxym = fSnxy/fNHits;
      nx2m = fSnx2/fNHits;
      ny2m = fSny2/fNHits;
      nxrm = fSnxr/fNHits;
      nyrm = fSnyr/fNHits;

      double d  = nx2m*ny2m-nxym*nxym;

      double xc = (nyrm*nx2m-nxrm*nxym)/d;
      double yc = (nyrm*nxym-nxrm*ny2m)/d;

      CofM.SetX(xc);
      CofM.SetY(yc);
//-----------------------------------------------------------------------------
// calculate seed chi2 - can this be optimized ?
//-----------------------------------------------------------------------------
      fChi2All  = 0;
      fChi2Perp = 0;

      for (int face=0; face<kNFaces; face++) {
        const HitData_t* hd = hitData[face];
        if (hd) {
          double dx = hd->fHit->pos().x()-xc;
          double dy = hd->fHit->pos().y()-yc;
//-----------------------------------------------------------------------------
// split into wire parallel and perpendicular components
//-----------------------------------------------------------------------------
          const XYZVectorF& wdir = hd->fHit->wdir();

          double dxy_dot_w = dx*wdir.x()+dy*wdir.y();
          double dxy_dot_n = dx*wdir.y()-dy*wdir.x();

          float  chi2_par  = (dxy_dot_w*dxy_dot_w)/(SigmaR2+hd->fSigW2);
          float  chi2_perp = (dxy_dot_n*dxy_dot_n)/SigmaR2;
          float  chi2      = chi2_par + chi2_perp;
          fChi2All        += chi2;
          fChi2Perp       += chi2_perp;
        }
      }
    }
//-----------------------------------------------------------------------------
// utility: calculate Com and chi2's, call for N>= 2 hit seeds
//-----------------------------------------------------------------------------
  void DeltaSeed::Chi2(double Xc, double Yc, double SigmaR2, double& Chi2All, double& Chi2Perp) {

    Chi2All  = 0;
    Chi2Perp = 0;

    for (int face=0; face<kNFaces; face++) {
      const HitData_t* hd = hitData[face];
      if (hd) {
        double dx = hd->fHit->pos().x()-Xc;
        double dy = hd->fHit->pos().y()-Yc;
//-----------------------------------------------------------------------------
// split into wire parallel and perpendicular components
//-----------------------------------------------------------------------------
        const XYZVectorF& wdir = hd->fHit->wdir();

        double dxy_dot_w    = dx*wdir.x()+dy*wdir.y();
        double dxy_dot_n    = dx*wdir.y()-dy*wdir.x();

        float  chi2_par     = (dxy_dot_w*dxy_dot_w)/(SigmaR2+hd->fSigW2);
        float  chi2_perp    = (dxy_dot_n*dxy_dot_n)/SigmaR2;
        float  chi2         = chi2_par + chi2_perp;
        Chi2All            += chi2;
        Chi2Perp           += chi2_perp;
      }
    }
  }
}
