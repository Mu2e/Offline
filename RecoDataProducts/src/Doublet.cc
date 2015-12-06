//
#include "RecoDataProducts/inc/Doublet.hh"


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
    fNStrawHits = 1;
    
    fTrkDir[0]  = trkdir;
    fTrkPos[0]  = trkpos;
    fHit   [0]  = hit;
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
