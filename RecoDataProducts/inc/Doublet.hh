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
#include "BTrk/KalmanTrack/KalRep.hh"
#include "KalmanTests/inc/TrkDef.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "BTrk/TrkBase/HelixParams.hh"

// C++
#include <vector>

using namespace std; 

namespace mu2e {

  class Doublet{
  public:
    enum { kMaxHits = 5 } ;

    int                 fIndex;		// doublet index in the list
    int                 fStationId;
    int                 fPanelId;
    CLHEP::Hep3Vector   fShDir;
    int                 fNStrawHits;
    int                 fStrawAmbig[kMaxHits];
    CLHEP::Hep3Vector   fTrkDir    [kMaxHits];
    CLHEP::Hep3Vector   fTrkPos    [kMaxHits];
    mu2e::TrkStrawHit*  fHit       [kMaxHits];
    int                 fHitIndex  [2]; // indices of the used pair of hits 
    double              fTrkDxDz;	// track dx/dz in the local panel frame
    double              fDxDz      [4];	// 4 combinations
    double              fChi2      [4];	// 4 combinations
    int                 fIBest;		// best combination
    int                 fINext;		// next-to-best combination
    int                 fOs;		// 0:opposite sign, +/-2:same sign
    double              fMcDoca    [kMaxHits]; // signed MC distance of closest approach
//-----------------------------------------------------------------------------
// constructors and such
//-----------------------------------------------------------------------------
    Doublet();

    Doublet(int index,
	    int station, int panel, 
	    CLHEP::Hep3Vector shdir, 
	    CLHEP::Hep3Vector trkdir,
	    CLHEP::Hep3Vector trkpos,
	    TrkStrawHit*      hit);
    
    ~Doublet();

    double Chi2Best() const { return fChi2[fIBest]; }
    
    void addStrawHit(CLHEP::Hep3Vector trkdir,
		     CLHEP::Hep3Vector trkpos,
		     TrkStrawHit*      hit);

    int    isSameSign() const { return (fIBest == 0) || (fIBest == 2) ; }

    double bestDxDz   () const { return fDxDz[fIBest];         }
    double trkDxDz    () const { return fTrkDxDz;              }
    double bestDxDzRes() const { return fTrkDxDz-fDxDz[fIBest];}
    
  };
  
  typedef std::vector<Doublet> DoubletCollection;
}

#endif
