//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////


#include "Offline/CalPatRec/inc/DeltaFinder_types.hh"

namespace mu2e {
  namespace DeltaFinderTypes {

//-----------------------------------------------------------------------------
    void Data_t::printHitData(HitData_t* Hd, const char* Option) {
      const mu2e::ComboHit* ch0 = &(*chcol)[0];

      const mu2e::ComboHit* ch = Hd->fHit;

      int index = ch-ch0;

      printf("index   sid  Stn:Pln:Pnl:Str seed delta  Time     TCorr     eDep    chi2min     sigw\n");
      printf("------------------------------------------\n");
      printf("%5i %5i ",
             index,ch->strawId().asUint16());

      printf(" %3i %3i %3i %3i",
             ch->strawId().station(),
             ch->strawId().plane(),
             ch->strawId().panel(),
             ch->strawId().straw());

      printf("%5i %5i",Hd->fSeedIndex,Hd->fDeltaIndex);
      printf("  %8.3f %8.3f %8.5f",ch->time(),ch->correctedTime(),ch->energyDep());
      printf(" %8.2f %8.2f\n",Hd->fChi2Min,Hd->fSigW);
    }

//-----------------------------------------------------------------------------
    void Data_t::printDeltaSeed(DeltaSeed* Seed, const char* Option) {
      printf("---------------------------------------------");
      printf("-------------------------------------------------------------------------------------\n");
      printf("index good type delta  SHID  SHID  SHID  SHID");
      printf("  chi2all/N    chi21   chi22 mintime  maxtime     X        Y         Z   nfwh nch nsh\n");
      printf("---------------------------------------------");
      printf("-------------------------------------------------------------------------------------\n");

      printf("%5i  %03i %4i %5i",Seed->Index(),Seed->fGood,Seed->fType,Seed->fDeltaIndex);
//-----------------------------------------------------------------------------
// print hit ID's in each face
//-----------------------------------------------------------------------------
      for (int face=0; face<kNFaces; face++) {
        const HitData_t* hd = Seed->hitData[face];
        if (hd == nullptr) printf(" %5i",-1);
        else {
          const ComboHit* hit = hd->fHit;
          printf(" %5i",hit->strawId().asUint16());
        }
      }

      printf(" %8.2f %7.2f %7.2f",Seed->Chi2AllDof(),Seed->fChi21,Seed->fChi22);
      printf("%8.1f %8.1f",Seed->MinHitTime(),Seed->MaxHitTime());
      printf(" %8.3f %8.3f %9.3f",Seed->CofM.x(),Seed->CofM.y(),Seed->CofM.z());
      printf("%4i",Seed->fNFacesWithHits);
      printf("%4i",Seed->NHits());
      printf("%4i",Seed->NStrawHits());
      printf("\n");
    }

//-----------------------------------------------------------------------------
    void Data_t::printDeltaCandidate(DeltaCandidate* Delta, const char* Option) {

      printf("------------------------------------------------------------------------------------------------------\n");
      printf("      i    nh  ns s1  s2     X       Y         Z          chi21   chi22   htmin   htmax   t0min   t0max     \n");
      printf("------------------------------------------------------------------------------------------------------\n");
      printf(":dc:%05i %3i",Delta->Index(),Delta->fNHits);
      printf(" %3i",Delta->n_seeds);
      printf(" %2i  %2i %7.2f %7.2f %9.2f",Delta->fFirstStation,Delta->fLastStation,
             Delta->CofM.x(),Delta->CofM.y(),Delta->CofM.z());
      printf("\n");
      printf("------------------------------------------------------------------------------------------------------\n");

      for (int is=Delta->fFirstStation;is<=Delta->fLastStation; is++) {
        DeltaSeed* ds = Delta->seed[is];
        if (ds != NULL) {

          int face0 = ds->SFace(0);
          int face1 = ds->SFace(1);

          const HitData_t* hd0 = ds->HitData(face0);
          const HitData_t* hd1 = (face1 >= 0) ? ds->HitData(face1) : nullptr;

          printf("          %3i  %3i    %3i:%03i",ds->fNHits,ds->fNHitsCE,is,ds->Index());
          printf(" %7.2f %7.2f %9.2f",ds->CofM.x(),ds->CofM.y(),ds->CofM.z());
          float chi22 = (hd1) ? hd1->fChi2Min : -1;
          printf(" %7.1f %7.1f",hd0->fChi2Min, chi22);
          printf(" %7.1f %7.1f",ds->MinHitTime(),ds->MaxHitTime());
          printf(" %7.1f %7.1f",Delta->fT0Min[is]  ,Delta->fT0Max[is]);

          printf("  (");
          for (int face=0; face<kNFaces; face++) {
            const HitData_t* hd = ds->HitData(face);
            if (hd == nullptr) printf(" %5i",-1);
            else {
              const ComboHit* hit = hd->fHit;
              printf(" %5i",hit->strawId().asUint16());
            }
            if (face != kNFaces-1) printf(",");
          }

          printf(")\n");
        }
      }
    }

//-----------------------------------------------------------------------------
    DeltaSeed::DeltaSeed() {
      printf("ERROR: DeltaSeed::DeltaSeed() should not be called\n");
    }

//-----------------------------------------------------------------------------
    DeltaSeed::DeltaSeed(int Index, int Station, int Face0, HitData_t* Hd0, int Face1, HitData_t* Hd1) {
      fIndex            = Index;
      fStation          = Station;
      fGood             =  1;
      fSFace[0]         = Face0;
      fSFace[1]         = Face1;

      for (int face=0; face<kNFaces; face++) {
        fFaceProcessed[face] = 0;
        hitData       [face] = NULL;
        fMcPart       [face] = NULL;
      }

      hitData[Face0]    = Hd0;

      fChi21            = Hd0->fChi2Min;

      if (Face1 >= 0) {
        hitData[Face1]    = Hd1;
        fType             =  10*Face0+Face1;
        fNHits            =  2;
        fNStrawHits       = Hd0->fHit->nStrawHits()+Hd1->fHit->nStrawHits();
        fNFacesWithHits   =  2;
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

      fDeltaIndex       = -1;
      fChi2All          = 9999.99;
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
    }

//-----------------------------------------------------------------------------
    DeltaCandidate::DeltaCandidate() {
      fIndex  = -1;
      for(int s=0; s<kNStations; ++s) {
        dxy    [s] = -1;
        seed   [s] = NULL;
        fT0Min [s] = -1.e10;
        fT0Max [s] =  1.e10;
      }
      fFirstStation = 999;
      fLastStation  =  -1;
      fMcPart       = NULL;
      fNHits        = 0;
      fNStrawHits   = 0;
      fNHitsCE      = 0;
      n_seeds       = 0;
    }

    DeltaCandidate::DeltaCandidate(int Index, DeltaSeed* Seed, int Station) {
      fIndex  = Index;
      for(int s=0; s<kNStations; ++s) {
        dxy    [s] = -1;
        seed   [s] = NULL;
        fT0Min [s] = -1.e10;
        fT0Max [s] =  1.e10;
      }
      fFirstStation = 999;
      fLastStation  =  -1;
      fMcPart       = NULL;
      fNHits        = 0;
      fNStrawHits   = 0;
      fNHitsCE      = 0;
      n_seeds       = 0;

      if (Seed) AddSeed(Seed,Station);
    }

//-----------------------------------------------------------------------------
    void DeltaCandidate::AddSeed(DeltaSeed* Seed, int Station) {
      seed[Station]         = Seed;
      n_seeds              += 1;

      if (fFirstStation > Station) fFirstStation = Station;
      if (fLastStation  < Station) fLastStation  = Station;
//-----------------------------------------------------------------------------
// ** FIXME recalculate the center of gravity - dont' need to be exact here
// seeds with NHits=1 don't have their center of gravity defined
//-----------------------------------------------------------------------------
      if (Seed->fNFacesWithHits > 1) {
        float x = (CofM.x()*fNHits+Seed->CofM.x()*Seed->NHits())/(fNHits+Seed->NHits());
        float y = (CofM.y()*fNHits+Seed->CofM.y()*Seed->NHits())/(fNHits+Seed->NHits());
        CofM.SetX(x);
        CofM.SetY(y);
      }
                                        // precalculate phi
      phi                   = CofM.phi();
      fNHits               += Seed->NHits();
      fNStrawHits          += Seed->NStrawHits();
//-----------------------------------------------------------------------------
// update t0min and t0max
// FIXME - need more reasonable limits
//-----------------------------------------------------------------------------
      float t0min = Seed->T0Min();
      float t0max = Seed->T0Max();
      float t0    = (t0min+t0max)/2;
      float dt    = 50-(t0max-t0min)/2;
      fT0Min[Station] = t0-dt;
      fT0Max[Station] = t0+dt;
//-----------------------------------------------------------------------------
// finally, set seed DeltaIndex
//-----------------------------------------------------------------------------
      Seed->fDeltaIndex = fIndex;
    }

//-----------------------------------------------------------------------------
    void DeltaCandidate::MergeDeltaCandidate(DeltaCandidate* Delta,int PrintErrorDiagnostics) {
      int is1 = Delta->FirstStation();
      int is2 = Delta->LastStation ();
      for (int is=is1; is<=is2; is++) {
        if (Delta->seed[is] == nullptr)                               continue;
        if ((seed[is] != nullptr) and PrintErrorDiagnostics) {
          printf("ERROR in DeltaCandidate::MergeDeltaCandidate: ");
          printf("merged DC also has a segment in station %i\n",is);
                                                                      continue;
        }
        seed[is] = Delta->seed[is];
        seed[is]->fDeltaIndex = fIndex;
//-----------------------------------------------------------------------------
// increment hit count only if a seed has been addded
//-----------------------------------------------------------------------------
        n_seeds              += 1;
        fNHits               += Delta->seed[is]->NHits();
        fNStrawHits          += Delta->seed[is]->NStrawHits();

        if (fFirstStation > is) fFirstStation = is;
        if (fLastStation  < is) fLastStation  = is;
      }

      float x = (CofM.x()*fNHits+Delta->CofM.x()*Delta->NHits())/(fNHits+Delta->NHits());
      float y = (CofM.y()*fNHits+Delta->CofM.y()*Delta->NHits())/(fNHits+Delta->NHits());
      CofM.SetX(x);
      CofM.SetY(y);
    }

//-----------------------------------------------------------------------------
    int findIntersection(const HitData_t* Hd1, const HitData_t* Hd2, Intersection_t* Result) {
      double x1, y1, x2, y2, nx1, ny1, nx2, ny2;

      const XYZVectorF& p1 = Hd1->fHit->pos();
      x1 =  p1.x();
      y1 =  p1.y();

      const XYZVectorF& p2 = Hd2->fHit->pos();
      x2 =  p2.x();
      y2 =  p2.y();

      const XYZVectorF& wdir1 = Hd1->fHit->wdir();
      nx1 = wdir1.x();
      ny1 = wdir1.y();

      const XYZVectorF& wdir2 = Hd2->fHit->wdir();
      nx2 = wdir2.x();
      ny2 = wdir2.y();

      double n1n2  = nx1*nx2+ny1*ny2;
      double r12n1 = (x1-x2)*nx1+(y1-y2)*ny1;
      double r12n2 = (x1-x2)*nx2+(y1-y2)*ny2;
//-----------------------------------------------------------------------------
// intersection point, in 2D two lines always intersect
// wd1, wd2 - distances to the hits, sign convention: delta=hit-intersection
//-----------------------------------------------------------------------------
      double t1   = (r12n2*n1n2-r12n1)/(1-n1n2*n1n2);
      Result->wd1 = -t1;

      Result->x = x1+nx1*t1;
      Result->y = y1+ny1*t1;
      Result->z = (p1.z()+p2.z())/2;

      double t2   = (r12n2-n1n2*r12n1)/(1-n1n2*n1n2);
      Result->wd2 = -t2;

      return 0;
    }
  }
}
