#include "CalPatRec/inc/Doublet.hh"


namespace mu2e {
  
  Doublet::Doublet(int index,
		   int station, int panel, 
		   CLHEP::Hep3Vector shdir, 
		   CLHEP::Hep3Vector trkdir,
		   CLHEP::Hep3Vector trkpos,
		   TrkStrawHit*      hit){
    fDoubletIndex = index;
    fStationId    = station; 
    fPanelId      = panel;
    fShDir        = shdir;
    fNstrawHits   = 1;
    
    fTrkDirCol[0] = trkdir;
    fTrkPosCol[0] = trkpos;
    fTrkshcol [0] = hit;
  }
  
  Doublet::~Doublet(){}

  void Doublet::addStrawHit(CLHEP::Hep3Vector trkdir,
			    CLHEP::Hep3Vector trkpos,
			    TrkStrawHit*      hit){
    if (fNstrawHits <5){
    //set the parameters of the new staw hit
    fTrkDirCol[fNstrawHits] = trkdir;
    fTrkPosCol[fNstrawHits] = trkpos;
    fTrkshcol [fNstrawHits] = hit;

    //now increment the size 
    ++fNstrawHits;
    } else {
      printf("[Doublet::addStrawhit] ERROR: trying to add more than 5 hits in the same panel!\n");
    }
  }
}
