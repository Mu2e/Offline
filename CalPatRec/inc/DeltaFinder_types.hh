#ifndef __CalPatRec_DeltaFinder_types_hh__
#define __CalPatRec_DeltaFinder_types_hh__

namespace art {
  class Event;
};

namespace fhicl {
  class ParameterSet;
};

#include "DataProducts/inc/StrawId.hh"
#include "RecoDataProducts/inc/StereoHit.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "TrackerGeom/inc/Tracker.hh"

namespace mu2e {
  class Panel;
  class SimParticle;
//-----------------------------------------------------------------------------
// delta-electron seed: structure within the station
// doesn't own anything, no need to delete any pinters
//-----------------------------------------------------------------------------
  namespace DeltaFinderTypes {
//-----------------------------------------------------------------------------
// intersection of the two hit wires
//-----------------------------------------------------------------------------
    struct Intersection_t {
      double     x;			// x-coordinate of the intersection point
      double     y;			// y-coordinate of the intersection point
      double     t1;			// distanCe from the center of the 1st wire
      double     t2;			// distance from the center of the 2nd wire
      double     wd1;                   // distance btw the 1st hit and the intersection
      double     wd2;		        // distance btw the 2nd hit and the intersection
    };

    struct DeltaCandidate;
    
    enum {
      kNStations      = StrawId::_nplanes/2,   // number of tracking stations
      kNFaces         = StrawId::_nfaces*2 ,   // N(faces) per station (4)
      kNPanelsPerFace = StrawId::_npanels/2    // = 3
    };
    
    struct HitData_t {
      const ComboHit*         fHit;
      // const StrawHitPosition* fPos;
      // const Straw*            fStraw;
      int                     fSeedNumber;
      int                     fNSecondHits;
      int                     fDeltaIndex;
      float                   fChi2Min;
      float                   fSigW;     // cached resolution along the wire
      float                   fRMid;
      float                   fDr;	// work variable

      HitData_t(const ComboHit* Hit, /*const StrawHitPosition* Pos, const Straw* aStraw,*/ float SigW) {
	fHit         = Hit; 
	// fPos         = Pos; 
	// fStraw       = aStraw; 
	fChi2Min     = 1.1e10; 
	fSigW        = SigW; 
	fSeedNumber  = -1; 
	fNSecondHits = -1;
	fDeltaIndex  = -1;
	fRMid        = sqrt(fHit->centerPos().Mag2());//FIXME! crosscheck  fStraw->getMidPoint().perp();
	fDr          = 1.1e10;
      }

      int Used() { return (fChi2Min < 1.e10) ; }
    };

    struct PanelZ_t {
      int                              fNHits  ; // guess, total number of ComboHits
      std::vector<HitData_t>           fHitData;
      const Panel*                     fPanel;      // backward pointer to the tracker panel
      double                           wx;          // direction cosines of the wires, assumed to be all the same
      double                           wy;
      double                           phi;         // phi angle of the wire
      double                           z;           // 
    }; 

//-----------------------------------------------------------------------------
// diagnostics structure
//-----------------------------------------------------------------------------
    struct McPart_t {
      const SimParticle*            fSim;        // type undefined here
      const DeltaCandidate*         fDelta;
      std::vector<const HitData_t*> fListOfHits;
      int                           fFirstStation;
      int                           fLastStation;
      int                           fID;	           // SimParticle::id().asInt()
      int                           fPdgID;
      int                           fNHitsCE;
      int                           fNHitsDelta;	   // number of hits associated with reconstructed delta electrons
      float                         fTime;
      float                         fHitTMin;          // min and max times of the straw hits
      float                         fHitTMax;
      float                         fStartMom;

      McPart_t(const SimParticle* Sim = NULL) { 
	fSim          = Sim; 
	fDelta        = NULL;
	fID           = -1;
	fPdgID        = 0;
	fNHitsDelta   = 0;
	fNHitsCE      = 0;
	fFirstStation = 999;
	fLastStation  = -1;
	fStartMom     = -1;
	fTime         = 1.e6; 		// at initialization, make it absurd
	fHitTMin      = 1.e6;
	fHitTMax      = -1.e6;
      }

      ~McPart_t() { fListOfHits.clear(); }

      int NHits() { return fListOfHits.size(); }

      float Momentum () { return fStartMom; }
      float Time     () { return fTime; }
      float HitDt    () { return fHitTMax-fHitTMin; }
    }; 

    class DeltaSeed {
    public:
      int                            fNumber;		// number within the station
      int                            fStation;		// station
      int                            fType;             // defines indices of the two faces used for preseed seach
					                // 0: used in recovery
      int                            fFaceProcessed[kNFaces];
      int                            fGood;             // <-killer number> if not to be used - what about 0 ?
      int                            fNFacesWithHits;
      int                            fNHitsTot;         // total number of hits
      int                            fNHitsCE;
      const HitData_t*               fHitData[2];       // stereo hit seeding the seed search
      // const StrawHitPosition*        fPos[2];
      float                          fChi21;             // chi2's of the two initial hits
      float                          fChi22;
      PanelZ_t*                      panelz   [kNFaces];
      std::vector<const HitData_t*>  hitlist  [kNFaces];
      std::vector<McPart_t*>         fMcPart  [kNFaces];
      CLHEP::Hep3Vector              CofM;
      float                          fMinTime;          // min and max times of the included hits
      float                          fMaxTime;
      float                          fMaxDriftTime;
      McPart_t*                      fPreSeedMcPart[2]; // McPart_t's corresponding to initial intersection 
      int                            fDeltaIndex;
      float                          fChi2All;

      DeltaSeed() {
	fType             =  0;
	fGood             =  1;
	fNHitsTot         =  0;
	fNHitsCE          =  0;
	fNFacesWithHits   =  0;
	fStation          = -1;
	fNumber           = -1;
	fPreSeedMcPart[0] = NULL;
	fPreSeedMcPart[1] = NULL;
	//	used              = false;
	fMinTime          = 1.e10;
	fMaxTime          = -1.;
	fMaxDriftTime     = -1.;
	for (int face=0; face<kNFaces; face++) {
	  fFaceProcessed[face] = 0;
	  panelz        [face] = NULL;
	}
	fDeltaIndex       = -1;
	fChi21            = -1;
	fChi22            = -1;
	fChi2All          = 1.e10;
      }
      //------------------------------------------------------------------------------
      // dont need a copy constructor
      //------------------------------------------------------------------------------
      ~DeltaSeed() {}

      float            Chi2N   ()         { return (fChi21+fChi22)/fNHitsTot; }
      float            Chi2Tot ()         { return (fChi21+fChi22); }
      float            Chi2All ()         { return fChi2All; }
      float            Chi2AllDof ()      { return fChi2All/fNHitsTot; }
      int              NHits   (int Face) { return hitlist[Face].size(); }
      const HitData_t* HitData (int Face, int I) { return hitlist[Face][I]; } // no boundary check !
      int              NHitsTot()         { return fNHitsTot; }
      int              MCTruth ()         { return (fPreSeedMcPart[0] != NULL) && (fPreSeedMcPart[0] == fPreSeedMcPart[1]) ; }
      bool             Used    ()         { return (fDeltaIndex >= 0); }
      //-----------------------------------------------------------------------------
      // drift time can't be < 0, so it is the minimal measured time which represents 
      // best the 'particle time'
      // keep in mind that fMaxTime < particle T0 < fMinTime
      //-----------------------------------------------------------------------------
      float            Time()     { return fMinTime; }
      float            T0Min()    { return fMaxTime-fMaxDriftTime;   }
      float            T0Max()    { return fMinTime; }
    };

    struct DeltaCandidate {
    public:
      int                   fNumber;
      int                   fFirstStation;
      int                   fLastStation;
      DeltaSeed*            seed   [kNStations];
      float                 dxy    [kNStations];   // used only for diagnostics
      float                 fT0Min [kNStations];
      float                 fT0Max [kNStations];
      CLHEP::Hep3Vector     CofM;
      float                 phi;
      int                   n_seeds;
      McPart_t*             fMcPart;
      int                   fNHits;
      int                   fNHitsMcP;	       // Nhits by the "best" particle"
      int                   fNHitsCE;

      DeltaCandidate() {
	fNumber = -1;
	for(int s=0; s<kNStations; ++s) {
	  dxy    [s] = -1;
	  seed   [s] = NULL;
	}
	fFirstStation = -1;
	fLastStation  = -1;
	fMcPart       = NULL;
	fNHits        = 0;
	fNHitsCE      = 0;
      }

      int        NSeeds               () { return n_seeds; }
      DeltaSeed* Seed            (int I) { return seed[I]; }
      bool       StationUsed     (int I) { return (seed[I] != NULL); }
      float      T0Min           (int I) { return fT0Min[I]; }
      float      T0Max           (int I) { return fT0Max[I]; }
    };
//-----------------------------------------------------------------------------
// data structure passed to the histogramming routine
//-----------------------------------------------------------------------------
    struct Data_t {
      const art::Event*             event;
      const Tracker*                tracker;
      std::string                   strawDigiMCCollectionTag;
      std::string                   ptrStepPointMCVectorCollectionTag;
      std::vector<DeltaSeed*>       seedHolder [kNStations];
      std::vector<DeltaCandidate>   deltaCandidateHolder;
      PanelZ_t                      oTracker[kNStations][kNFaces][kNPanelsPerFace];
      int                           stationUsed[kNStations];
      int                           nseeds;
      int                           nseeds_per_station[kNStations];
      const ComboHitCollection*     chcol;
      const StrawHitFlagCollection* shfcol;
      const TimeClusterCollection*  tpeakcol;
      int                           debugLevel;	     // printout level
    };

//-----------------------------------------------------------------------------
// finally, utility functions
//-----------------------------------------------------------------------------
    int findIntersection(const HitData_t* Hd1, const HitData_t* Hd2, Intersection_t* Result);
  }
}
#endif
