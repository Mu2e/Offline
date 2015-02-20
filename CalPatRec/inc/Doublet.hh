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
    Doublet(int station, int panel, 
	    CLHEP::Hep3Vector shdir, 
	    CLHEP::Hep3Vector trkdir,
	    CLHEP::Hep3Vector trkpos,
	    TrkStrawHit *hit):
      fStationId(station), fPanelId(panel),
      fShDir(shdir){
      fTrkDirCol.push_back(trkdir);
      fTrkPosCol.push_back(trkpos);
      fTrkshcol.push_back(hit);
    }
    
    ~Doublet(){}
    
    int fStationId;
    int fPanelId;
    CLHEP::Hep3Vector fShDir;
    std::vector<CLHEP::Hep3Vector> fTrkDirCol;
    std::vector<CLHEP::Hep3Vector> fTrkPosCol;

    std::vector< mu2e::TrkStrawHit*> fTrkshcol;
    

  };
  
  typedef std::vector<Doublet> DoubletCollection;
}

#endif
