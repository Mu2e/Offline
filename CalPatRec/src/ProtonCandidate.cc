//
#include "Offline/CalPatRec/inc/ProtonCandidate.hh"

namespace mu2e {
//-----------------------------------------------------------------------------
  ProtonCandidate::ProtonCandidate(int Index) :
    fSt(0), fSz(0), fSt2(0), fStz(0), fSz2(0)
  {
    fIndex             = Index;
    fFirstStation      = 999;
    fLastStation       =  -1;
    fNStationsWithHits = 0;
    fNHits             = 0;
    fNStrawHits        = 0;
    fSumEDep           = 0;

    fMcPart            = nullptr;
    fNHitsCE           = 0;
  }

//-----------------------------------------------------------------------------
// first added seed has at least one stereo, so the COG calculation COG is safe
//-----------------------------------------------------------------------------
  void ProtonCandidate::addSeed(DeltaSeed* Seed) {
    int station = Seed->Station();
    float z     = DeltaFinderTypes::stationZ[station];

    if (fFirstStation > station) fFirstStation = station;
    if (fLastStation  < station) fLastStation  = station;
//-----------------------------------------------------------------------------
// ** FIXME recalculate the center of gravity - dont' need to be exact here
// seeds with SFace(1) < 0 don't have their center of gravity defined - their hits
// have been picked up individually, the intersection doesn't matter
//-----------------------------------------------------------------------------
    // int nh = Seed->NHits();

    for (int face=0; face<kNFaces; face++) {
      HitData_t* hd = Seed->HitData(face);
      if (hd == nullptr)                                              continue;
//-----------------------------------------------------------------------------
// now we neeeed to figure whether this hit has already been added
//-----------------------------------------------------------------------------
      int found = 0;
      int nh = fHitData[station][face].size();
      for (int ih=0; ih<nh; ih++) {
        HitData_t* hd2 = fHitData[station][face][ih];
        if (hd == hd2) {
          found = 1;
          break;
        }
      }

      if (found == 0) {
        if (fNHitsStation[station] == 0) fNStationsWithHits++;
        fNHitsStation[station] += 1;

        fHitData[station][face].push_back(hd);
        fNHits      += 1;
        fNStrawHits += hd->fHit->nStrawHits();
        fSumEDep    += hd->fHit->energyDep();

        float t      = hd->fCorrTime;

        fSt         += t;
        fSt2        += t*t;
        fStz        += t*z;
        fSz         += z;
        fSz2        += z*z;

      }
    }

    // CofM.SetX(xc);
    // CofM.SetY(yc);
//-----------------------------------------------------------------------------
// recalculate the timing parameterization
//-----------------------------------------------------------------------------
    double tm, zm, tzm, t2m, z2m;

    tm  = fSt /fNHits;
    zm  = fSz /fNHits;
    tzm = fStz/fNHits;
    t2m = fSt2/fNHits;
    z2m = fSz2/fNHits;
                                        // 'combo-hit'-based way, FIXME
    if (fNStationsWithHits > 1) {
      fDtDz  = (tzm-tm*zm)/(z2m-zm*zm);
      fT0    = tm-fDtDz*zm;
      fSigT0 = sqrt((t2m-tm*tm)/(fNHits-1+1.e-12));
    }
    else {
      fDtDz  = 0;
      fT0    = tm;
      fSigT0 = 0;
    }
  }

//-----------------------------------------------------------------------------
// first added seed has at least one stereo, so the COG calculation COG is safe
// by default, UpdateTime=1
//-----------------------------------------------------------------------------
  void ProtonCandidate::addHit(int Station, HitData_t* Hd, int UpdateTime) {

//-----------------------------------------------------------------------------
// the hit has been already included, do nothing
//-----------------------------------------------------------------------------
    int face = Hd->fZFace;
    for (auto & hd: fHitData[Station][face]) { if (hd == Hd) return; }

    if (fNHitsStation[Station] == 0) fNStationsWithHits++;
    fNHitsStation[Station] += 1;

    if (fFirstStation > Station) fFirstStation = Station;
    if (fLastStation  < Station) fLastStation  = Station;

    if (fPanelID[Station][face] == -1) fPanelID[Station][face] = Hd->panelID();
//-----------------------------------------------------------------------------
// ** FIXME recalculate the center of gravity - dont' need to be exact here
// seeds with SFace(1) < 0 don't have their center of gravity defined - their hits
// have been picked up individually, the intersection doesn't matter
//-----------------------------------------------------------------------------
    fNHits      += 1;
    fNStrawHits += Hd->fHit->nStrawHits();
    fSumEDep    += Hd->fHit->energyDep();
//-----------------------------------------------------------------------------
// time: the sums need to be updated always
//-----------------------------------------------------------------------------
    double t = Hd->fCorrTime;
    double z = DeltaFinderTypes::stationZ[Station];
    fSt     += t;
    fSt2    += t*t;
    fStz    += t*z;
    fSz     += z;
    fSz2    += z*z;
//-----------------------------------------------------------------------------
// by default, update the time .. but may want to skip this step
//-----------------------------------------------------------------------------
    if (UpdateTime) {
      double tm, zm, tzm, t2m, z2m;

      tm  = fSt /fNHits;
      zm  = fSz /fNHits;
      tzm = fStz/fNHits;
      t2m = fSt2/fNHits;
      z2m = fSz2/fNHits;
                                        // 'combo-hit'-based way, FIXME
      if (fNStationsWithHits > 1) {
        fDtDz  = (tzm-tm*zm)/(z2m-zm*zm);
        fT0    = tm-fDtDz*zm;
        fSigT0 = sqrt((t2m-tm*tm)/(fNHits-0.9999));
      }
      else {
        fDtDz  = 0;
        fT0    = tm;
        fSigT0 = 0;
      }
    }
  }

//-----------------------------------------------------------------------------
  float ProtonCandidate::Phi(int Station) {
    printf("ProtonCandidate::%s not impemented yet\n",__func__);
    return -100;
  }

}
