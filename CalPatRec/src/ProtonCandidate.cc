//
#include "Offline/CalPatRec/inc/ProtonCandidate.hh"

namespace mu2e {
//-----------------------------------------------------------------------------
// fIndex is set just once upon construction
//-----------------------------------------------------------------------------
  ProtonCandidate::ProtonCandidate(int Index) {
    fIndex = Index;
    init();
  }

//-----------------------------------------------------------------------------
  void ProtonCandidate::init() {
    fMask              = 0;
    fFirstStation      = 999;
    fLastStation       =  -1;
    fNStationsWithHits = 0;
    fNHitsTot          = 0;
    fNStrawHitsTot     = 0;
    fSumEDep           = 0.;
    fTMid              = 0.;

    fT0                = 0.f;
    fDtDz              = 0.f;
    fSigT0             = 0.f;

    fMcPart            = nullptr;
    fNHitsMcP          = 0;
    fNHitsCE           = 0;
    fTimeIndex         = 0;

    for (int is=0; is<kNStations; is++) {
      fNHitsStation[is] = 0;
      fMinHitTime  [is] = 0;
      fMaxHitTime  [is] = 0;
      fSumX        [is] = 0.;
      fSumY        [is] = 0.;
      fPhi         [is] = 0.;
      fNHitsStation[is] = 0.;
      for (int face=0; face<kNFaces; face++) {
        fHitData[is][face].clear();
        fPanelID[is][face] = -1;
      }
    }

    fSt  = 0.;
    fSz  = 0.;
    fSt2 = 0.;
    fStz = 0.;
    fSz2 = 0.;
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

        if (fHitData[station][face].size() == 0) {
//-----------------------------------------------------------------------------
// adding first hit
//-----------------------------------------------------------------------------
          fPanelID[station][face] = hd->panelID();
        }

        if (fPanelID[station][face] == hd->panelID()) {
//-----------------------------------------------------------------------------
// adding not first hit
//-----------------------------------------------------------------------------
          fHitData[station][face].push_back(hd);
          fNHitsTot      += 1;
          fNStrawHitsTot += hd->fHit->nStrawHits();
          fSumEDep       += hd->fHit->energyDep()*hd->fHit->nStrawHits();

          fSumX[station] += hd->fX;
          fSumY[station] += hd->fY;

          float t      = hd->fCorrTime;

          fSt         += t;
          fSt2        += t*t;
          fStz        += t*z;
          fSz         += z;
          fSz2        += z*z;
        }
        else {
//-----------------------------------------------------------------------------
// dont' know exactly what to do, for the moment - print diagnostics
//-----------------------------------------------------------------------------
          printf("ProtonCandidate::%s ERROR: trying to add hit in a wrong panel\n",__func__);
        }
      }
    }
//-----------------------------------------------------------------------------
// knowing phi may be helpful
//-----------------------------------------------------------------------------
    fPhi[station] = atan2(fSumY[station],fSumX[station]);
//-----------------------------------------------------------------------------
// recalculate the timing parameterization
//-----------------------------------------------------------------------------
    updateTime();
//-----------------------------------------------------------------------------
// update hit indices
//-----------------------------------------------------------------------------
    Seed->setProtonIndex(fIndex);
    for (int face=0; face<kNFaces; face++) {
      if (Seed->HitData(face)) Seed->HitData(face)->setProtonIndex(fIndex);
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

    fHitData[Station][face].push_back(Hd);

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
    fNHitsTot      += 1;
    fNStrawHitsTot += Hd->fHit->nStrawHits();
    fSumEDep       += Hd->fHit->energyDep()*Hd->fHit->nStrawHits();
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
    if (UpdateTime) updateTime();
    Hd->setProtonIndex(fIndex);
  }

//-----------------------------------------------------------------------------
// assume hit is there
//-----------------------------------------------------------------------------
  void ProtonCandidate::removeHit(int Station, HitData_t* Hd, int UpdateTime) {
//-----------------------------------------------------------------------------
// the hit has been already included, do nothing
//-----------------------------------------------------------------------------
    int face = Hd->fZFace;
    int nhf = fHitData[Station][face].size();
    for (int ih=0; ih<nhf; ih++) {
      HitData_t* hd = fHitData[Station][face][ih];
      if (hd == Hd) {
        fNHitsStation[Station] -= 1;
        if (fNHitsStation[Station] == 0) fNStationsWithHits--;

        if (fNHitsStation[Station] == 0) {
//-----------------------------------------------------------------------------
// may need to redefined the first or the last station
//-----------------------------------------------------------------------------
          if (fFirstStation == Station) {
            for (int is=Station+1; is<=fLastStation; is++) {
              if (fNHitsStation[is] > 0) {
                fFirstStation = is;
                break;
              }
            }
          }
          else if (Station == fLastStation) {
            for (int is=Station-1; is>=fFirstStation; is--) {
              if (fNHitsStation[is] > 0) {
                fLastStation = is;
                break;
              }
            }
          }
        }
//-----------------------------------------------------------------------------
// ** FIXME recalculate the center of gravity - dont' need to be exact here
// seeds with SFace(1) < 0 don't have their center of gravity defined - their hits
// have been picked up individually, the intersection doesn't matter
//-----------------------------------------------------------------------------
        fNHitsTot      -= 1;
        int nsh         = Hd->fHit->nStrawHits();
        fNStrawHitsTot -= nsh;
        fSumEDep       -= Hd->fHit->energyDep()*nsh;
//-----------------------------------------------------------------------------
// time: the sums need to be updated always
//-----------------------------------------------------------------------------
        double t = Hd->fCorrTime;
        double z = DeltaFinderTypes::stationZ[Station];
        fSt     -= t;
        fSt2    -= t*t;
        fStz    -= t*z;
        fSz     -= z;
        fSz2    -= z*z;
//-----------------------------------------------------------------------------
// finally, remove hit from the hit list. That may require moving other hits
// If face becomes empty, flag that
//-----------------------------------------------------------------------------
        for (int ih2=ih+1; ih2<nhf; ih2++) {
          fHitData[Station][face][ih2-1] = fHitData[Station][face][ih2];
        }
        fHitData[Station][face].pop_back();

        if (nhf == 1) fPanelID[Station][face] = -1;
      }
    }
//-----------------------------------------------------------------------------
// by default, update the time .. but may want to skip this step
//-----------------------------------------------------------------------------
    if (UpdateTime) updateTime();
  }

//-----------------------------------------------------------------------------
  void ProtonCandidate::predictPhi(int Station, PhiPrediction_t* Prediction) {

    Prediction->fPhi     = -100;
    Prediction->fErr     =  0;
    Prediction->fPanelID = -1;
//-----------------------------------------------------------------------------
// cant predict phi with one station
//-----------------------------------------------------------------------------
    if (Station < fFirstStation) {
//-----------------------------------------------------------------------------
// first pass... - moving upstream
//-----------------------------------------------------------------------------
// can't predict anything if the gap is too large
//-----------------------------------------------------------------------------
      if (fFirstStation-Station > 1)                                    return;
//-----------------------------------------------------------------------------
// prediction right into the next station, make sure thre are hits in the two
// previous ones
//-----------------------------------------------------------------------------
      if (fNStationsWithHits > 1) {
//-----------------------------------------------------------------------------
// there are hits in the two last consecutive stations
//-----------------------------------------------------------------------------
        float phi1       = fPhi[fLastStation ];
        float phi2       = fPhi[fFirstStation];
        float dphi       = (phi2-phi1);
        if (dphi < M_PI) dphi += 2*M_PI;
        if (dphi > M_PI) dphi -= 2*M_PI;

        dphi             = dphi/(fLastStation-fFirstStation);
        Prediction->fPhi = phi2+dphi;
        Prediction->fErr = 0.5; // fmax(0.2,fabs(dphi));
      }
      else {
//-----------------------------------------------------------------------------
// just one station
//-----------------------------------------------------------------------------
        float phi = fPhi[fFirstStation  ];
        Prediction->fPhi = phi+0.15;
        Prediction->fErr = 0.65;
      }
    }
    else {
//-----------------------------------------------------------------------------
// picking hits up in between the first and the last stations
// don't know what to do yet
//-----------------------------------------------------------------------------
      if (fNHitsStation[Station] > 0) {
//-----------------------------------------------------------------------------
// one of the stations with hits
//-----------------------------------------------------------------------------
        Prediction->fPhi = phi(Station);
        Prediction->fErr = 0.2;
      }
      else {
//-----------------------------------------------------------------------------
// if a proton candidate has enough hits, one could use a parabola **FIXME**
//-----------------------------------------------------------------------------
        float phi1       = fPhi[fLastStation ];
        float phi2       = fPhi[fFirstStation];
        float dphi       = (phi2-phi1)/(fLastStation-fFirstStation);
        Prediction->fPhi = phi2-dphi*(Station-fFirstStation);
        Prediction->fErr = 0.5; // fmax(0.2,fabs(dphi));
      }
    }
  }

//-----------------------------------------------------------------------------
  void ProtonCandidate::updateTime() {

    double tm, zm, tzm, t2m, z2m;

    tm  = fSt /fNHitsTot;
    zm  = fSz /fNHitsTot;
    tzm = fStz/fNHitsTot;
    t2m = fSt2/fNHitsTot;
    z2m = fSz2/fNHitsTot;
                                        // 'combo-hit'-based way, FIXME ? or just OK ?
    if (fNStationsWithHits > 1) {
      fDtDz  = (tzm-tm*zm)/(z2m-zm*zm);
      fT0    = tm-fDtDz*zm;
      fSigT0 = sqrt((t2m-tm*tm)/(fNHitsTot-1));
    }
    else {
      fDtDz  = 0;
      fT0    = tm;
      fSigT0 = 0;
    }

    int station(9);

    if (fFirstStation <= 9) {
      if (fLastStation <= 9) station = fLastStation;
    }
    else {
      station = fFirstStation;
    }

    float z = DeltaFinderTypes::stationZ[station];
    fTMid   = fT0+fDtDz*z;
  }
}
