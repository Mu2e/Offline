///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"

#include "Offline/CalPatRec/inc/DeltaFinderAlg.hh"

namespace mu2e {
//-----------------------------------------------------------------------------
// start building proton candidates from high-eDep delta seeds
//-----------------------------------------------------------------------------
  int  DeltaFinderAlg::createProtonCandidates() {
    if(_doTiming > 0) _watch->SetTime(__func__);
    int kMaxProtonGap(1);                                // inefficiency gap
    PhiPrediction_t pred;                                // for phi predictions

    for (int station=kNStations-1; station>=0; station--) {

      int nseeds = _data->nProtonSeeds(station);
      for (int iseed=0; iseed<nseeds; iseed++) {
        DeltaSeed* seed = _data->protonSeed(station,iseed);
//-----------------------------------------------------------------------------
// skip segments considered to be duplicates of others
// fGood<-2000 : proton segment candidates
//-----------------------------------------------------------------------------
        if ((seed->fGood < 0) and (seed->fGood > -2000))              continue;
//-----------------------------------------------------------------------------
// some of the 'proton' seeds could be already flagged as delta's,
// skip them and do not try to count them as protons
//-----------------------------------------------------------------------------
        if (seed->deltaIndex() >= 0)                                  continue;
        float t         = seed->TMean();
        int   found{0};
//-----------------------------------------------------------------------------
// check if 'seed' is consistent in time with any of existing proton candidates
// the number of proton candidates could increase during the processing
// of the previous seed
//-----------------------------------------------------------------------------
        int npc = _data->nProtonCandidates();
        for (int ipc=0; ipc<npc; ipc++) {
          ProtonCandidate* pc = _data->protonCandidate(ipc);
          if (pc->fFirstStation-station >  kMaxProtonGap)             continue;
          float tpc = pc->t0(station);
          if (fabs(tpc-t) > 30)                                       continue;
//-----------------------------------------------------------------------------
// the proton candidate and the seed are close in time, see if can compare phi
// if predicted phi = -100, prediction couldn't be made
//-----------------------------------------------------------------------------
          int consistent(1);
          pc->predictPhi(station,&pred);
          if (pred.fPhi > -99) {
//-----------------------------------------------------------------------------
// phi prediction exists, check
//-----------------------------------------------------------------------------
            consistent = 0;

            float dphi = pred.fPhi-seed->CofM.phi();
            if (dphi < -M_PI) dphi += 2*M_PI;
            if (dphi >  M_PI) dphi -= 2*M_PI;

            if (fabs(dphi) < pred.fErr)  consistent = 1;
          }
//-----------------------------------------------------------------------------
// add hits of the seed to the proton candidate. Remember that there could
// seeds with the overlapping hits
// the same segment can be associated with several proton candidates -
// shooting for efficiency, don't need to be exclusive
//-----------------------------------------------------------------------------
          if (consistent) {
            pc->addSeed(seed);
            found +=1;
          }
        }
//-----------------------------------------------------------------------------
// if the seed has been associated with an existing proton candidate,
// proceed with the next seed
//-----------------------------------------------------------------------------
        if (found == 0) {
//-----------------------------------------------------------------------------
// at this point, the number of proton candidates is incremented by one
//-----------------------------------------------------------------------------
          ProtonCandidate* pc = _data->newProtonCandidate();
          pc->init();
          pc->addSeed(seed);
        }
      }
    }

    if(_doTiming > 0) _watch->StopTime(__func__);
    return 0;
  }

//-----------------------------------------------------------------------------
// merge proton candidates with overlapping hit patterns
//-----------------------------------------------------------------------------
  int  DeltaFinderAlg::resolveProtonOverlaps(std::vector<ProtonCandidate*>* Pc) {
    if(_doTiming > 1) _watch->SetTime(__func__);
//-----------------------------------------------------------------------------
// assume no more than 100 overlapping hits per face
//-----------------------------------------------------------------------------
    constexpr int max_overlap = 100;
    HitData_t* ohit[kNStations][kNFaces][max_overlap];
    int        novr[kNStations][kNFaces];

    int nprot = Pc->size();

    for (int i1=0; i1<nprot-1; i1++) {
      ProtonCandidate* p1 = Pc->at(i1);
      float t1 = p1->tMid();
      for (int i2=i1+1; i2<nprot; i2++) {
        ProtonCandidate* p2 = Pc->at(i2);
        float t2 = p2->tMid();
//-----------------------------------------------------------------------------
// by construction, p2->time >= p1->time()
//-----------------------------------------------------------------------------
        float dt = t2-t1;

        int os1 = std::max(p1->fFirstStation,p2->fFirstStation);
        int os2 = std::min(p1->fLastStation ,p2->fLastStation );
        int ds  = os2-os1;
        if (ds < 0) ds=0;

        if (dt > 60 + 5*ds)                                           break;
//-----------------------------------------------------------------------------
// the two candidates are close enough in time, count overlapping hits
//-----------------------------------------------------------------------------
        int noverlap = 0;

        for (int is=os1; is<=os2; is++) {
          int nh1 = p1->nHitsStation(is);
          if (nh1 == 0)                                               continue;
          int nh2 = p2->nHitsStation(is);
          if (nh2 == 0)                                               continue;
          for (int face=0; face<kNFaces; face++) {
            novr[is][face] = 0;
            int nhf1 = p1->nHits(is,face);
            int nhf2 = p2->nHits(is,face);
            if ((nhf1 == 0) or (nhf2 == 0))                           continue;
            for (int ih1=0; ih1<nhf1; ih1++) {
              HitData_t* h1 = p1->hitData(is,face,ih1);
              for (int ih2=0; ih2<nhf2; ih2++) {
                HitData_t* h2 = p2->hitData(is,face,ih2);
                if (h1 == h2) {
                  const int novr_fz = novr[is][face];
                  if(novr_fz >= max_overlap) { // FIXME: Should this be an error?
                    // throw std::runtime_error("Too many overlapping hits!");
                    continue;
                  }
                  ohit[is][face][novr_fz] = h1;
                  novr[is][face] += 1;
                  noverlap       += 1;
                  break;
                }
              }
            }
          }
        }

        int nh1tot = p1->nHitsTot();
        int nh2tot = p2->nHitsTot();

        if      (noverlap >= 2) {
          if (nh2tot >= nh1tot) {
//-----------------------------------------------------------------------------
// merge p1 into p2 by removing overlapping hits from p1
//-----------------------------------------------------------------------------
            for (int is=os1; is<=os2; is++) {
              for (int face=0; face<kNFaces; face++) {
                int nhf = novr[is][face];
                for (int ih=0; ih<nhf; ih++) {
                  HitData_t* h = ohit[is][face][ih];
                  p1->removeHit(is,h,0);
                  h->setProtonIndex(p2->index());
                }
              }
            }
            p1->updateTime();
            if (_debugLevel > 0) {
              printf("* DeltaFinderAlg::%s:001: overlapping hits hits proton segments %3i:%3i, delete them from %3i\n",
                     __func__,p1->index(),p2->index(),p1->index());
            }
            break;
          }
          else {
//-----------------------------------------------------------------------------
// merge p2 into p1
//-----------------------------------------------------------------------------
            for (int is=os1; is<=os2; is++) {
              for (int face=0; face<kNFaces; face++) {
                int nhf = novr[is][face];
                for (int ih=0; ih<nhf; ih++) {
                  HitData_t* h = ohit[is][face][ih];
                  p2->removeHit(is,h,0);
                  h->setProtonIndex(p1->index());
                }
              }
            }
            p2->updateTime();
            if (_debugLevel > 0) {
              printf("* DeltaFinderAlg::%s:001: overlapping hits hits proton segments %3i:%3i, delete them from %3i\n",
                     __func__,p1->index(),p2->index(),p2->index());
            }
          }
        }
      }
    }

    if(_doTiming > 1) _watch->StopTime(__func__);
    return 0;
  }

//-----------------------------------------------------------------------------
  int  DeltaFinderAlg::mergeNonOverlappingCandidates(std::vector<ProtonCandidate*>* Pc) {
    if(_doTiming > 1) _watch->SetTime(__func__);

    int nprot = Pc->size();
    PhiPrediction_t pred; // for phi predictions

    for (int i1=0; i1<nprot-1; i1++) {
      ProtonCandidate* p1 = Pc->at(i1);
      if (p1->nStationsWithHits() < 2)                                continue;
      float t1 = p1->tMid();
      for (int i2=i1+1; i2<nprot; i2++) {
        ProtonCandidate* p2 = Pc->at(i2);
        if (p2->nStationsWithHits() < 2)                              continue;
        float t2 = p2->tMid();
//-----------------------------------------------------------------------------
// by construction, p2->time >= p1->time()
//-----------------------------------------------------------------------------
        float dt = t2-t1;
        int os1  = std::max(p1->fFirstStation,p2->fFirstStation);
        int os2  = std::min(p1->fLastStation ,p2->fLastStation );
        int gap  = os1-os2;
//------------------------------------------------------------------------------
// require the two candidates to be close enough in time
//-----------------------------------------------------------------------------
        if (gap > 0) {
//-----------------------------------------------------------------------------
// segments do not overlap, assume proton propagation speed 5 ns/station
//-----------------------------------------------------------------------------
          float t1, t1_pred;

          if (p1->fFirstStation > p2->fLastStation) {
            t1 = p1->t0(p1->fFirstStation);
            float t2 = p2->t0(p2->fLastStation);
            t1_pred = t2 + 5*(p1->fFirstStation-p2->fLastStation);
          }
          else {
            t1 = p1->t0(p1->fLastStation);
            float t2 = p2->t0(p2->fFirstStation);
            t1_pred = t2 - 5*(p2->fFirstStation-p1->fLastStation);
          }

          if (t1_pred-t1 > 40)                                          break;
//-----------------------------------------------------------------------------
// check gap (in units of stations) between the candidates
// 'pl' - long, 'ps' - short
//-----------------------------------------------------------------------------
          ProtonCandidate  *pl(p1), *ps(p2);

          int l1 = p1->fLastStation-p1->fFirstStation+1;
          int l2 = p2->fLastStation-p2->fFirstStation+1;

          int lmax = std::max(l1,l2);
          if (l2 > l1) {
            pl = p2;
            ps = p1;
          }
//-----------------------------------------------------------------------------
// if gap is too large, can't reliably extrapolate and predict phi, so do nothing
//-----------------------------------------------------------------------------
          if      (lmax < gap)                                        continue;
//-----------------------------------------------------------------------------
// gap is not too large, 'p' is the longest, use it to predict phi
// figure the station to predict it to
//-----------------------------------------------------------------------------
          int s0 = ps->fLastStation;
          if (pl->fLastStation < ps->fFirstStation) s0 = ps->fFirstStation;
          pl->predictPhi(s0,&pred);
          float phi_s = ps->phi(s0);
          float dphi  = pred.fPhi-phi_s;
          if      (dphi >  M_PI) dphi -= 2*M_PI;
          else if (dphi < -M_PI) dphi += 2*M_PI;

          if (fabs(dphi) > 0.3)                                       continue;
          if (_debugLevel > 0) {
            printf("* close non-overlapping in Z segments segments: p1, p2: %5i %5i \n",
                   p1->index(),p2->index());
          }
        }
        else {
//-----------------------------------------------------------------------------
// segments do overlap
//-----------------------------------------------------------------------------
          if (dt > 30.f)                                                 break;
//-----------------------------------------------------------------------------
// segments are close in time, check phi in station=os1 (in the overlap region)
//-----------------------------------------------------------------------------
          float phi1 = p1->phi(os1);
          float phi2 = p2->phi(os1);
          if (fabs(phi1-phi2) > 0.3f)                                 continue;
//-----------------------------------------------------------------------------
// the two segments are close in time- see how many false positives are there..
// not many, can use...
//-----------------------------------------------------------------------------
          if (_debugLevel > 0) {
            printf("* close     overlapping in Z segments segments: p1, p2: %5i %5i \n",
                   p1->index(),p2->index());
          }
        }
      }
    }

    if(_doTiming > 1) _watch->StopTime(__func__);
    return 0;
  }

//-----------------------------------------------------------------------------
  int  DeltaFinderAlg::mergeProtonCandidates() {
    if(_doTiming > 1) _watch->SetTime(__func__);
//-----------------------------------------------------------------------------
// to speed up the execution, start from time-ordering the list
//-----------------------------------------------------------------------------
    int nprot = _data->nProtonCandidates();

    std::vector<ProtonCandidate*> pc;
    pc.reserve(nprot);
    for (int i=0; i<nprot; i++) pc.push_back(_data->protonCandidate(i));

    std::sort(pc.begin(),pc.end(),
              [](const ProtonCandidate* a, const ProtonCandidate* b) { return a->fTMid < b->fTMid; });
//-----------------------------------------------------------------------------
// for debugging purposes, set index corresponding to time ordering
//-----------------------------------------------------------------------------
    for (int i=0; i<nprot; i++) {
      pc[i]->fTimeIndex = i;
    }
//-----------------------------------------------------------------------------
// step one: resolve overlaps: if two proton candidates have the same segment,
// assigne it to the one with more hits
// usually, it happens with candiates which are 2-station long,
// and it is the right thing to do
//-----------------------------------------------------------------------------
    resolveProtonOverlaps(&pc);
//-----------------------------------------------------------------------------
// step two: merge (supposedly, non-overlapping) parts of the same proton trajectory
// consider only candidates which have segments in at least two stations
//-----------------------------------------------------------------------------
    if (_mergePC) mergeNonOverlappingCandidates(&pc);

    if(_doTiming > 1) _watch->StopTime(__func__);
    return 0;
  }

//-----------------------------------------------------------------------------
  int  DeltaFinderAlg::findProtons() {
    if(_doTiming > 0) _watch->SetTime(__func__);
//-----------------------------------------------------------------------------
// fill data structure with "proton" hits only (E > 3 keV)
    prepareProtonHits();
//-----------------------------------------------------------------------------
// loop over all stations and find proton time clusters
// start from seeds
//-----------------------------------------------------------------------------
    createProtonCandidates();
//-----------------------------------------------------------------------------
// merge proton candidates with overlapping hit lists
//-----------------------------------------------------------------------------
    mergeProtonCandidates();
//-----------------------------------------------------------------------------
// last step - check stations without found segments and try to pick up missing hits
//-----------------------------------------------------------------------------
    if (_pickupProtonHits) recoverMissingProtonHits();

    if(_doTiming > 0) _watch->StopTime(__func__);
    return 0;
  }


//-----------------------------------------------------------------------------
// don't expect protons to have just one hit per station in many stations
// proton candidate is alowed to have hits only in one panel per face
//-----------------------------------------------------------------------------
  int  DeltaFinderAlg::recoverMissingProtonHits() {

    int npc = _data->nProtonCandidates();
    PhiPrediction_t pred; // for phi predictions
    for (int ipc=0; ipc<npc; ipc++) {
      ProtonCandidate* pc = _data->protonCandidate(ipc);
//-----------------------------------------------------------------------------
// don't extend candidates made out of one segment - but there is no such
// start from the first station to define limits
//-----------------------------------------------------------------------------
      int s1 = std::max(pc->fFirstStation-1,0);
      int s2 = std::min(pc->fLastStation+1,(int)kNStations-1);
//-----------------------------------------------------------------------------
// check inside "holes"
//-----------------------------------------------------------------------------
      for (int station=s2; station>=s1; station--) {
//-----------------------------------------------------------------------------
// predict proton time and try to predict phi, very roughly
//-----------------------------------------------------------------------------
        float proton_time = pc->t0(station);

        pc->predictPhi(station,&pred);

        for (int face=0; face<kNFaces; face++) {
          FaceZ_t* fz = _data->faceData(station,face);
//-----------------------------------------------------------------------------
// consider only faces where the candidate so far has no hits
//-----------------------------------------------------------------------------
          if (pc->nHits(station,face) > 0)                            continue;

          int nph = fz->nProtonHits();
          for (int iph=0; iph<nph; iph++) {
            HitData_t* hd = fz->protonHitData(iph);
            if (fabs(hd->fCorrTime-proton_time) > _maxDriftTime)      continue;
            if (pred.fPhi > -99) {
//-----------------------------------------------------------------------------
// an idea of phi prediction exists, check the hit phi
//-----------------------------------------------------------------------------
              float dphi = pred.fPhi-hd->Phi();
              if (dphi < -M_PI) dphi += 2*M_PI;
              if (dphi >  M_PI) dphi -= 2*M_PI;

              if (fabs(dphi) > pred.fErr)                             continue;
            }
            int panel_id = pc->fPanelID[station][face];
            if ((panel_id == -1) or (hd->panelID() == panel_id)) {
//-----------------------------------------------------------------------------
// ** FIXME** here we add the first hit possible - this is a room for a mistake - SURE
// add hit to the list of proton hits, does the hit need to be marked ? - YES
// what do we do in case of confusion ? - NOTHING, keep the last one
// hit knows about its Z-face.. a hit could be a part of more than one time cluster,
// so need to check explicitly
//-----------------------------------------------------------------------------
              pc->addHit(station,hd);
              if (_debugLevel > 0) {
                printf("DeltaFinderAlg::%s: add combo hit sid=%5i station:face: %2i:%i to pc:%3i\n",
                       __func__,hd->fHit->strawId().asUint16(),station,face,pc->index());
              }
            }
          }
        }
      }
    }

    return 0;
  }

//-----------------------------------------------------------------------------
  int  DeltaFinderAlg::prepareProtonHits() {

    for (int is=0; is<kNStations; is++) {
      for (int face=0; face<kNFaces; face++) {
        FaceZ_t* fz = _data->faceData(is,face);
        int nh = fz->nHits();
        for (int ih=0; ih<nh; ih++) {
          HitData_t* hd = fz->hitData(ih);
//-----------------------------------------------------------------------------
// require hit to have EDep consistent with 95+% of the protons and not to
// be a part of a reconstructed delta electron candidate
// this should significantly reduce the number of candidate hits to consider
//-----------------------------------------------------------------------------
          if ((hd->fHit->energyDep() > _minProtonHitEDep) and (hd->fDeltaIndex < 0)) {

            fz->fProtonHitData.push_back(hd);

            // int nph = fz->nProtonHits();
            // int tbin = int (hd->fHit->time()/_timeBin) ;
            // if (fz->fPFirst[tbin] < 0) fz->fPFirst[tbin] = nph;
            // fz->fPLast[tbin] = nph;
          }
        }
      }
    }

    return 0;
  }
}
