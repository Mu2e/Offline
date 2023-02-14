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

    for (int station=kNStations-1; station>=0; station--) {

      int nseeds = _data->nProtonSeeds(station);
      for (int iseed=0; iseed<nseeds; iseed++) {
        DeltaSeed* seed = _data->protonSeed(station,iseed);
        float t         = seed->TMean();
        float rho       = seed->Rho();
        float seed_nx   = seed->Xc()/rho;
        float seed_ny   = seed->Yc()/rho;
        int   found{0};
//-----------------------------------------------------------------------------
// check if 'seed' is consistent in time with any of existing proton candidates
// the number of proton candidates could increase during the processing
// of the previous seed
//-----------------------------------------------------------------------------
        int npc = _data->nProtonCandidates();
        for (int ipc=0; ipc<npc; ipc++) {
          ProtonCandidate* pc = _data->protonCandidate(ipc);
          float tpc = pc->T0(station);
          if (fabs(tpc-t) > 20)                                       continue;
//-----------------------------------------------------------------------------
// the proton candidate and the seed are close in time, see if can compare phi
// if predicted phi = -100, prediction couldn't be made
//-----------------------------------------------------------------------------
          float phi_pc = pc->Phi(station);
          if (phi_pc < -99)                                           continue;
//-----------------------------------------------------------------------------
// prediction exists, check
//-----------------------------------------------------------------------------
          float nx = cos(phi_pc);
          float ny = sin(phi_pc);

          float n1n2 = nx*seed_nx+ny*seed_ny;
          if (n1n2 < 0)                                               continue;
//-----------------------------------------------------------------------------
// add hits of the seed to the proton candidate. Remember that there could
// seeds with the overlapping hits
//-----------------------------------------------------------------------------
          pc->addSeed(seed);
          found = 1;
                                                                      break;
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
          pc->addSeed(seed);
        }
      }
    }

    return 0;
  }

//-----------------------------------------------------------------------------
  int  DeltaFinderAlg::findProtons() {
//-----------------------------------------------------------------------------
// fill data structure with "proton" hits only (E > 3 keV)
    prepareProtonHits();
//-----------------------------------------------------------------------------
// loop over all stations and find proton time clusters
// start from seeds
//-----------------------------------------------------------------------------
    createProtonCandidates();
//-----------------------------------------------------------------------------
// last step - check stations without found segments and try to pick up missing hits
//-----------------------------------------------------------------------------
    recoverMissingProtonHits();

    return 0;
  }


//-----------------------------------------------------------------------------
// don't expect protons to have just one hit per station in many stations
// proton candidate is alowed to have hits only in one panel per face
//-----------------------------------------------------------------------------
  int  DeltaFinderAlg::recoverMissingProtonHits() {

    int npc = _data->nProtonCandidates();
    for (int ipc=0; ipc<npc; ipc++) {
      ProtonCandidate* pc = _data->protonCandidate(ipc);
//-----------------------------------------------------------------------------
// don't extend candidates made out of one segment - but there is no such
// start from the first station to define limits
//-----------------------------------------------------------------------------
      int s1 = pc->fFirstStation;
      int s2 = pc->fLastStation;
//-----------------------------------------------------------------------------
// check inside "holes"
//-----------------------------------------------------------------------------
      for (int station=s2; station>=s1; station--) {
//-----------------------------------------------------------------------------
// predict proton time and try to predict phi, very roughly
//-----------------------------------------------------------------------------
        float proton_time = pc->T0(station);
        float proton_phi  = pc->Phi(station);

        for (int face=0; face<kNFaces; face++) {
          FaceZ_t* fz = _data->faceData(station,face);

          int nph = fz->nProtonHits();
          for (int iph=0; iph<nph; iph++) {
            HitData_t* hd = fz->protonHitData(iph);
            if (fabs(hd->fCorrTime-proton_time) > _maxDriftTime)     continue;
            if (proton_phi > -99) {
//-----------------------------------------------------------------------------
// an idea of phi prediction exists, check the hit phi
//-----------------------------------------------------------------------------
              float nx = cos(proton_phi);
              float ny = sin(proton_phi);

              float n1n2 = nx*hd->fX+ny*hd->fY;
              if (n1n2 < 0)                                          continue;
            }
            int panel_id = pc->fPanelID[station][face];
            if ((panel_id == -1) or (hd->panelID() == panel_id)) {
//-----------------------------------------------------------------------------
// ** FIXME** here we add the first hit possible - this is a room for a mistake
// add hit to the list of proton hits, does the hit need to be marked ?
// what do we do in case of confusion ?
// hit knows about its Z-face.. a hit could be a part of more than one time cluster,
// so need to check explicitly
//-----------------------------------------------------------------------------
              pc->addHit(station,hd);
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
// this shoudl significantly reduce the number of candidate hits to consider
//-----------------------------------------------------------------------------
          if ((hd->fHit->energyDep() > _minProtonHitEDep) and (hd->fDeltaIndex < 0)) {

            int nph = fz->nProtonHits();
            fz->fProtonHitData.push_back(hd);

            int tbin = int (hd->fHit->time()/_timeBin) ;
            if (fz->fPFirst[tbin] < 0) fz->fPFirst[tbin] = nph;
            fz->fPLast[tbin] = nph;
          }
        }
      }
    }

    return 0;
  }
}
