//
/*
#include "Offline/BTrkData/inc/Doublet.hh"


namespace mu2e {

//-----------------------------------------------------------------------------
// make sure neither number makes sense upon initialization
//-----------------------------------------------------------------------------
  Doublet::Doublet() {
    fIndex      = -1;
    fStationId  = -1;
    fPanelId    = -1;
    fShDir.set(0,0,0);
    fNStrawHits = -1.;
    for (int i=0; i<kMaxNHits; i++) {
      fStrawAmbig[i] = -999;
      fHit       [i] = NULL;
      fMcDoca    [i] = -999.;
      fTrkDir    [i].set(0.,0.,0);
      fTrkPos    [i].set(0.,0.,0);
    }

    fHitIndex[0] = -1;
    fHitIndex[1] = -1;
    fTrkDxDz     = -999.;

    for (int i=0; i<kMaxNComb; i++) {
      fDxDz[i] = -999.;
      fChi2[i] = -999.;
    }

    fIBest = -1;
    fINext = -1;
    fOs    = -999;
  }


//-----------------------------------------------------------------------------
  Doublet::Doublet(int               Index,
                   int               Station,
                   int               Panel,
                   CLHEP::Hep3Vector Shdir,
                   CLHEP::Hep3Vector Trkdir,
                   CLHEP::Hep3Vector Trkpos,
                   TrkStrawHit*      Hit) {

    fIndex       = Index;
    fStationId   = Station;
    fPanelId     = Panel;
    fShDir       = Shdir;
    fNStrawHits  = 1;

    for (int i=0; i<kMaxNHits; i++) {
      fStrawAmbig[i] = -999;
      fTrkDir    [i].set(0.,0.,0);
      fTrkPos    [i].set(0.,0.,0);
      fHit       [i] = NULL;
      fMcDoca    [i] = -999.;
    }

    fTrkDir[0]   = Trkdir;
    fTrkPos[0]   = Trkpos;
    fHit   [0]   = Hit;
                                        // make sure these are undefined
    fHitIndex[0] = -1;
    fHitIndex[1] = -1;
    fIBest       = -1;
    fINext       = -1;
    fOs          = -1;
  }

//-----------------------------------------------------------------------------
  void Doublet::addStrawHit(CLHEP::Hep3Vector trkdir,
                            CLHEP::Hep3Vector trkpos,
                            TrkStrawHit*      hit   ) {

    if (fNStrawHits < 0) {
      printf("[Doublet::addStrawHit] ERROR: fNStrawHits = %i kMaxNHits = %i!\n",fNStrawHits,kMaxNHits);
    }
    else if (fNStrawHits < kMaxNHits) {
                                        // set the parameters of the new staw hit
      fTrkDir[fNStrawHits] = trkdir;
      fTrkPos[fNStrawHits] = trkpos;
      fHit   [fNStrawHits] = hit;
                                        // increment the size
      ++fNStrawHits;
    }
    // else {
    //   printf("[Doublet::addStrawHit] ERROR: trying to add more than %i hits in the same panel!\n",kMaxNHits);
    // }
  }
}
*/
