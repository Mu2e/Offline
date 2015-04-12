//
#include "KalmanTests/inc/Doublet.hh"


namespace mu2e {
  
  Doublet::Doublet() {
  }


  Doublet::Doublet(int index,
		   int station, int panel, 
		   CLHEP::Hep3Vector shdir, 
		   CLHEP::Hep3Vector trkdir,
		   CLHEP::Hep3Vector trkpos,
		   TrkStrawHit*      hit){
    fIndex      = index;
    fStationId  = station; 
    fPanelId    = panel;
    fShDir      = shdir;
    fNstrawHits = 1;
    
    fTrkDir[0]  = trkdir;
    fTrkPos[0]  = trkpos;
    fHit   [0]  = hit;
  }
  
  Doublet::~Doublet(){}

  void Doublet::addStrawHit(CLHEP::Hep3Vector trkdir,
			    CLHEP::Hep3Vector trkpos,
			    TrkStrawHit*      hit){
    if (fNstrawHits < 5) {
					// set the parameters of the new staw hit
    fTrkDir[fNstrawHits] = trkdir;
    fTrkPos[fNstrawHits] = trkpos;
    fHit   [fNstrawHits] = hit;
					// increment the size 
    ++fNstrawHits;
    } else {
      printf("[Doublet::addStrawhit] ERROR: trying to add more than 5 hits in the same panel!\n");
    }
  }
}
