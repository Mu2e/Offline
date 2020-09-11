///////////////////////////////////////////////////////////////////////////////
// utilities for the Module to perform BaBar Kalman fit
// 2015 - 02 - 17 G. Pezzullo created class for housing the straw hits doublets
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __CalPatRec_Doublet_hh__
#define __CalPatRec_Doublet_hh__

#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"

#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrk/TrkBase/HelixParams.hh"

// C++
#include <vector>

using namespace std; 

namespace mu2e {

  class Doublet{
  public:
    enum { kMaxNHits = 9,
	   kMaxNComb = 4
    } ;

    int                 fIndex;		// doublet index in the list
    int                 fStationId;
    int                 fPanelId;
    CLHEP::Hep3Vector   fShDir;
    int                 fNStrawHits;
    int                 fOldAmbig  [kMaxNHits];
    int                 fStrawAmbig[kMaxNHits];
    CLHEP::Hep3Vector   fTrkDir    [kMaxNHits];
    CLHEP::Hep3Vector   fTrkPos    [kMaxNHits];
    mu2e::TrkStrawHit*  fHit       [kMaxNHits];
    double              fMcDoca    [kMaxNHits]; // signed MC distance of closest approach
    int                 fHitIndex  [2];	        // indices of the used pair of hits 
    double              fTrkDxDz;	        // track dx/dz in the local panel frame
    double              fDxDz      [kMaxNComb]; // combinations
    double              fChi2Slope [kMaxNComb]; // 
    double              fChi2Coord [kMaxNComb]; // 
    double              fChi2      [kMaxNComb]; // 
    int                 fIBest;		        // best combination
    int                 fINext;		        // next best combination
    int                 fOs;		        // 0:opposite sign, +/-2:same sign
//-----------------------------------------------------------------------------
// constructors and such
//-----------------------------------------------------------------------------
    Doublet();
    Doublet(const Doublet& R);                 // unfortunately, need a copy constructor...

    Doublet(int index,
	    int station, int panel, 
	    CLHEP::Hep3Vector shdir, 
	    CLHEP::Hep3Vector trkdir,
	    CLHEP::Hep3Vector trkpos,
	    TrkStrawHit*      hit);
    
    ~Doublet();

    double chi2Best     () const { return fChi2     [fIBest]; }
    double chi2SlopeBest() const { return fChi2Slope[fIBest]; }
    double chi2CoordBest() const { return fChi2Coord[fIBest]; }
    
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
