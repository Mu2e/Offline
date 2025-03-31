//
#include "Offline/CalPatRec/inc/DeltaCandidate.hh"

namespace mu2e {
//-----------------------------------------------------------------------------
    DeltaCandidate::DeltaCandidate() :
      fNx(0.f), fNy(0.f), fSumEDep(0.f),
      fSnx2(0.), fSnxy(0.), fSny2(0.),fSnxr(0.),fSnyr(0.),
      fSt(0.), fSz(0.), fSt2(0.), fStz(0.), fSz2(0.),
      fT0(0.), fDtDz(0.), fSigT0(0.)
    {
      fIndex  = -1;
      fMask         = 0;
      fFirstStation = 999;
      fLastStation  =  -1;
      for(int s=0; s<kNStations; ++s) {
        fSeed   [s] = nullptr;
      }
      fNHits        = 0;
      fNStrawHits   = 0;
      fNSeeds       = 0;
    }

//-----------------------------------------------------------------------------
    DeltaCandidate::DeltaCandidate(int Index, DeltaSeed* Seed) :
      DeltaCandidate()
    {
      fIndex  = Index;
      if (Seed) AddSeed(Seed);
    }

//-----------------------------------------------------------------------------
// remove seed. Saved DeltaCandidate always has two seeds or more
//-----------------------------------------------------------------------------
  void DeltaCandidate::removeSeed(const DeltaCandidate* Delta, int Station) {
    *this = *Delta;
    DeltaSeed* seed = fSeed[Station];
    if (seed == nullptr)                                                return;
    fSeed[Station] = nullptr;
//-----------------------------------------------------------------------------
// general case: the 'seed' needs to be removed
//-----------------------------------------------------------------------------
    fNSeeds       -= 1;
    fNHits        -= seed->nHits();
    fNStrawHits   -= seed->nStrawHits();
//-----------------------------------------------------------------------------
// update the first and the last station numbers
//-----------------------------------------------------------------------------
    if (Station == fFirstStation) {
      while (fSeed[++fFirstStation] == nullptr) {}
    }

    if (Station == fLastStation) {
      while (fSeed[--fLastStation] == nullptr) {}
    }
//-----------------------------------------------------------------------------
// update coordinates
//-----------------------------------------------------------------------------
    fSumEDep    -= seed->SumEDep();

    fSnx2       -= seed->fSnx2;
    fSnxy       -= seed->fSnxy;
    fSny2       -= seed->fSny2;
    fSnxr       -= seed->fSnxr;
    fSnyr       -= seed->fSnyr;

    double d     = fSnx2*fSny2-fSnxy*fSnxy;
    double xc    = (fSnyr*fSnx2-fSnxr*fSnxy)/d;
    double yc    = (fSnyr*fSnxy-fSnxr*fSny2)/d;

    CofM.SetX(xc);
    CofM.SetY(yc);

    double rho = sqrt(xc*xc+yc*yc);
    fNx        = xc/rho;
    fNy        = yc/rho;
//-----------------------------------------------------------------------------
// update time
//-----------------------------------------------------------------------------
    double t = seed->TMean();
    double z = DeltaFinderTypes::stationZ[Station];
    fSt     -= t;
    fSt2    -= t*t;
    fStz    -= t*z;
    fSz     -= z;
    fSz2    -= z*z;
                                        // and update the time...
    double tm, zm, tzm, t2m, z2m;

    tm       = fSt /fNSeeds;
    zm       = fSz /fNSeeds;
    tzm      = fStz/fNSeeds;
    t2m      = fSt2/fNSeeds;
    z2m      = fSz2/fNSeeds;
                                        // 'combo-hit'-based way, FIXME
    if (fNSeeds > 1) {
      fDtDz  = (tzm-tm*zm)/(z2m-zm*zm);
      fT0    = tm-fDtDz*zm;
      fSigT0 = sqrt((t2m-tm*tm)/(fNSeeds-0.9999));
    }
    else {
      fDtDz  = 0;
      fT0    = tm;
      fSigT0 = 0;
    }
//-----------------------------------------------------------------------------
// do not update the seed
//-----------------------------------------------------------------------------
  }

//-----------------------------------------------------------------------------
// first added seed has at least one stereo, so the COG calculation COG is safe
//-----------------------------------------------------------------------------
  void DeltaCandidate::AddSeed(DeltaSeed* Seed) {
    int station    = Seed->Station();

    fSeed[station] = Seed;
    fNSeeds       += 1;

    if (fFirstStation > station) fFirstStation = station;
    if (fLastStation  < station) fLastStation  = station;
//-----------------------------------------------------------------------------
// ** FIXME recalculate the center of gravity - dont' need to be exact here
// seeds with SFace(1) < 0 don't have their center of gravity defined - their hits
// have been picked up individually, the intersection doesn't matter
//-----------------------------------------------------------------------------
    fNHits      += Seed->nHits();
    fNStrawHits += Seed->nStrawHits();
    fSumEDep    += Seed->SumEDep();

    fSnx2       += Seed->fSnx2;
    fSnxy       += Seed->fSnxy;
    fSny2       += Seed->fSny2;
    fSnxr       += Seed->fSnxr;
    fSnyr       += Seed->fSnyr;

    double d     = fSnx2*fSny2-fSnxy*fSnxy;
    double xc    = (fSnyr*fSnx2-fSnxr*fSnxy)/d;
    double yc    = (fSnyr*fSnxy-fSnxr*fSny2)/d;

    CofM.SetX(xc);
    CofM.SetY(yc);

    double rho = sqrt(xc*xc+yc*yc);
    fNx        = xc/rho;
    fNy        = yc/rho;
//-----------------------------------------------------------------------------
// time
//-----------------------------------------------------------------------------
    double t = Seed->TMean();
    double z = DeltaFinderTypes::stationZ[station];
    fSt     += t;
    fSt2    += t*t;
    fStz    += t*z;
    fSz     += z;
    fSz2    += z*z;
                                        // and update the time...
    double tm, zm, tzm, t2m, z2m;

    tm  = fSt /fNSeeds;
    zm  = fSz /fNSeeds;
    tzm = fStz/fNSeeds;
    t2m = fSt2/fNSeeds;
    z2m = fSz2/fNSeeds;
                                        // 'combo-hit'-based way, FIXME
    if (fNSeeds > 1) {
      fDtDz  = (tzm-tm*zm)/(z2m-zm*zm);
      fT0    = tm-fDtDz*zm;
      fSigT0 = sqrt((t2m-tm*tm)/(fNSeeds-0.9999));
    }
    else {
      fDtDz  = 0;
      fT0    = tm;
      fSigT0 = 0;
    }

  }

//-----------------------------------------------------------------------------
  void DeltaCandidate::MergeDeltaCandidate(DeltaCandidate* Delta,int PrintErrorDiagnostics) {
    int is1 = Delta->FirstStation();
    int is2 = Delta->LastStation ();
    for (int is=is1; is<=is2; is++) {
      if (Delta->fSeed[is] == nullptr)                                 continue;
      if ((fSeed[is] != nullptr) and PrintErrorDiagnostics) {
        printf("ERROR in DeltaCandidate::MergeDeltaCandidate: ");
        printf("merged DC also has a segment in station %i,",is);
        printf(" LEAVE THE EXISTING ONE\n");
                                                                      continue;
      }
      fSeed[is]              = Delta->fSeed[is];
      fSeed[is]->fDeltaIndex = fIndex;
//-----------------------------------------------------------------------------
// increment hit count only if a seed has been addded
//-----------------------------------------------------------------------------
      fNSeeds               += 1;
      fNHits                += fSeed[is]->nHits();
      fNStrawHits           += fSeed[is]->nStrawHits();
      fSumEDep              += fSeed[is]->SumEDep();

      if (fFirstStation > is) fFirstStation = is;
      if (fLastStation  < is) fLastStation  = is;
    }
//-----------------------------------------------------------------------------
// update XY sums
//-----------------------------------------------------------------------------
    fSnx2          += Delta->fSnx2;
    fSnxy          += Delta->fSnxy;
    fSny2          += Delta->fSny2;
    fSnxr          += Delta->fSnxr;
    fSnyr          += Delta->fSnyr;
//-----------------------------------------------------------------------------
// and recalculate the center-of-mass XY position
//-----------------------------------------------------------------------------
    double d        = fSnx2*fSny2-fSnxy*fSnxy;
    double xc       = (fSnyr*fSnx2-fSnxr*fSnxy)/d;
    double yc       = (fSnyr*fSnxy-fSnxr*fSny2)/d;

    CofM.SetX(xc);
    CofM.SetY(yc);

    double rho      = sqrt(xc*xc+yc*yc);
    fNx             = xc/rho;
    fNy             = yc/rho;
//-----------------------------------------------------------------------------
// time
//-----------------------------------------------------------------------------
    fSt     += Delta->fSt;
    fSt2    += Delta->fSt2;
    fStz    += Delta->fStz;
    fSz     += Delta->fSz;
    fSz2    += Delta->fSz2;

    double tm, zm, tzm, t2m, z2m;

    tm  = fSt /fNSeeds;
    zm  = fSz /fNSeeds;
    tzm = fStz/fNSeeds;
    t2m = fSt2/fNSeeds;
    z2m = fSz2/fNSeeds;
                                        // 'combo-hit'-based way, FIXME
    fDtDz  = (tzm-tm*zm)/(z2m-zm*zm);
    fT0    = tm-fDtDz*zm;
    fSigT0 = sqrt((t2m-tm*tm)/(fNSeeds-0.9999));
  }

//-----------------------------------------------------------------------------
  void DeltaCandidate::markHitsAsUsed() {
    for (int is=fFirstStation; is<=fLastStation; is++) {
      DeltaSeed* s = fSeed[is];
      if (s == nullptr)                                              continue;
      for (int face=0; face<kNFaces; face++) {
        HitData_t* hd = s->HitData(face);
        if (hd) hd->fDeltaIndex = fIndex;
      }
    }
  }
}
