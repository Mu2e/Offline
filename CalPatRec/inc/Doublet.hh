///////////////////////////////////////////////////////////////////////////////
// utilities for the Module to perform BaBar Kalman fit
// 2015 - 02 - 17 G. Pezzullo created class for housing the straw hits doublets
//
// $Id:  $
// $Author:  $
// $Date: $
///////////////////////////////////////////////////////////////////////////////
#ifndef __CalPatRec_Doublet_hh__
#define __CalPatRec_Doublet_hh__

#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"

#include "KalmanTests/inc/TrkDef.hh"
#include "KalmanTrack/KalRep.hh"
#include "KalmanTests/inc/TrkDef.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "TrkBase/HelixParams.hh"

// C++
#include <vector>

using namespace std; 

namespace mu2e {

  class Doublet{
  public:
    Doublet(int index,
	    int station, int panel, 
	    CLHEP::Hep3Vector shdir, 
	    CLHEP::Hep3Vector trkdir,
	    CLHEP::Hep3Vector trkpos,
	    TrkStrawHit*      hit);
    
    ~Doublet();
    
    void addStrawHit(CLHEP::Hep3Vector trkdir,
		     CLHEP::Hep3Vector trkpos,
		     TrkStrawHit*      hit);
    
    int                 fDoubletIndex;
    int                 fStationId;
    int                 fPanelId;
    CLHEP::Hep3Vector   fShDir;
    int                 fNstrawHits;
    int                 fStrawAmbig[5];
    CLHEP::Hep3Vector   fTrkDirCol [5];
    CLHEP::Hep3Vector   fTrkPosCol [5];
    mu2e::TrkStrawHit*  fTrkshcol  [5];

  };
  
  typedef std::vector<Doublet> DoubletCollection;
}

#endif
