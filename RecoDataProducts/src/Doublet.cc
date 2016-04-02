//
#include "RecoDataProducts/inc/Doublet.hh"


namespace mu2e {
  
  Doublet::Doublet() {
  }


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
  
  Doublet::~Doublet(){}

  void Doublet::addStrawHit(CLHEP::Hep3Vector trkdir,
			    CLHEP::Hep3Vector trkpos,
			    TrkStrawHit*      hit   ) {

    if (fNStrawHits < kMaxHits) {
					// set the parameters of the new staw hit
      fTrkDir[fNStrawHits] = trkdir;
      fTrkPos[fNStrawHits] = trkpos;
      fHit   [fNStrawHits] = hit;
					// increment the size 
      ++fNStrawHits;
    } 
    else {
      printf("[Doublet::addStrawHit] ERROR: trying to add more than %i hits in the same panel!\n",kMaxHits);
    }
  }
}
