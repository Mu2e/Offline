//
#include "BTrkData/inc/Doublet.hh"


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
// 
  Doublet::Doublet(const Doublet& R) {
    fIndex      = R.fIndex;
    fStationId  = R.fStationId;
    fPanelId    = R.fPanelId;
    fShDir      = R.fShDir;
    fNStrawHits = R.fNStrawHits;

    for (int i=0; i<kMaxNHits; i++) {
      fStrawAmbig[i] = R.fStrawAmbig[i];
      fHit       [i] = R.fHit[i];
      fMcDoca    [i] = R.fMcDoca[i];
      fTrkDir    [i] = R.fTrkDir[i];
      fTrkPos    [i] = R.fTrkPos[i];
    }
    
    fHitIndex[0] = R.fHitIndex[0];
    fHitIndex[1] = R.fHitIndex[1];
    fTrkDxDz     = R.fTrkDxDz;

    for (int i=0; i<kMaxNComb; i++) {
      fDxDz[i] = R.fDxDz[i];
      fChi2[i] = R.fChi2[i];
    }

    fIBest = R.fIBest;
    fINext = R.fINext;
    fOs    = R.fOs;
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
  Doublet::~Doublet(){}

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
