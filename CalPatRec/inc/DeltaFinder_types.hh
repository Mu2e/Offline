#ifndef __CalPatRec_DeltaFinder_types_hh__
#define __CalPatRec_DeltaFinder_types_hh__

namespace art {
  class Event;
};

namespace fhicl {
  class ParameterSet;
};

class TH1F;
class TH2F;

#include "CalPatRec/inc/LsqSums2.hh"
#include "RecoDataProducts/inc/StereoHit.hh"

namespace mu2e {
  class StrawHitPosition;
  class Panel;
  class SimParticle;
//-----------------------------------------------------------------------------
// delta-electron seed: structure within the station
// doesn't own anything, no need to delete any pinters
//-----------------------------------------------------------------------------
  namespace DeltaFinderTypes {
    struct DeltaCandidate;
    
    enum {
      kNStations      = 20,
      kNFaces         =  4,
      kNPanelsPerFace =  3
    };
    
    struct HitData_t {
      const StrawHit*         fHit;
      const StrawHitPosition* fPos;
      const Straw*            fStraw;
      int                     fSeedNumber;
      int                     fNSecondHits;
      float                   fChi2Min;
      float                   fSigW;     // cached resolution along the wire

      HitData_t(const StrawHit* Hit, const StrawHitPosition* Pos, const Straw* aStraw, float SigW) {
	fHit = Hit; fPos = Pos; fStraw = aStraw; fChi2Min = 1.1e10; fSigW = SigW; fSeedNumber = -1; fNSecondHits = -1;
      }

      int Used() { return (fChi2Min < 1.e10) ; }
    };

    struct PanelZ_t {
      int                              fNHits  [2]; // guess, total number of hits per layer
      std::vector<HitData_t>           fHitData[2];
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
      const SimParticle*           fSim;        // type undefined here
      const DeltaCandidate*        fDelta;
      std::vector<const StrawHit*> fListOfHits;
      int                          fFirstStation;
      int                          fLastStation;
      int                          fID;	          // SimParticle::id().asInt()
      int                          fPdgID;
      int                          fNHitsDelta;	// number of hits associated with reconstructed delta electrons
      float                        fTime;
      float                        fStartMom;

      McPart_t(const SimParticle* Sim = NULL) { 
	fSim          = Sim; 
	fDelta        = NULL;
	fID           = -1;
	fPdgID        = 0;
	fNHitsDelta   = 0;
	fFirstStation = 999;
	fLastStation  = -1;
	fStartMom     = -1;
	fTime         = 1.e6; 		// at initialization, make it absurd
      }

      ~McPart_t() { fListOfHits.clear(); }

      int NHits() { return fListOfHits.size(); }

      float Momentum() { return fStartMom; }
      float Time    () { return fTime; }
    }; 
//
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
      const HitData_t*               fHitData[2];       // stereo hit seeding the seed search
      const StrawHitPosition*        fPos[2];
      float                          chi2dof;           // for two initial hits
      PanelZ_t*                      panelz   [kNFaces];
      std::vector<const HitData_t*>  hitlist  [kNFaces];
      std::vector<McPart_t*>         fMcPart  [kNFaces];
      CLHEP::Hep3Vector              CofM;
      float                          fMinTime; // min and max times of the included hits
      float                          fMaxTime;
      McPart_t*                      fPreSeedMcPart[2]; // McPart_t's corresponding to initial intersection 
      bool                           used;   // denotes if seed has been used in connectseeds

      DeltaSeed() {
	fType             =  0;
	fGood             =  1;
	fNHitsTot         =  0;
	fNFacesWithHits   =  0;
	fStation          = -1;
	fNumber           = -1;
	fPreSeedMcPart[0] = NULL;
	fPreSeedMcPart[1] = NULL;
	used              = false;
	fMinTime          = 1.e10;
	fMaxTime          = -1.;
	for (int face=0; face<kNFaces; face++) {
	  fFaceProcessed[face] = 0;
	  panelz        [face] = NULL;
	}
      }
      //------------------------------------------------------------------------------
      // dont need a copy constructor
      //------------------------------------------------------------------------------
      ~DeltaSeed() {}

      float            Chi2    ()         { return chi2dof; }
      int              NHits   (int Face) { return hitlist[Face].size(); }
      const HitData_t* HitData (int Face, int I) { return hitlist[Face][I]; } // no boundary check !
      int              NHitsTot()         { return fNHitsTot; }
      int              MCTruth ()         { return (fPreSeedMcPart[0] != NULL) && (fPreSeedMcPart[0] == fPreSeedMcPart[1]) ; }
      //-----------------------------------------------------------------------------
      // drift time can't be < 0, so it is the minimal measured time which represents 
      // best the 'particle time'
      //-----------------------------------------------------------------------------
      float            Time()     { return fMinTime; }

    };

    struct DeltaCandidate {
    public:
      int                   fNumber;
      DeltaSeed*            seed   [kNStations];
      bool                  st_used[kNStations];
      float                 dxy    [kNStations];
      CLHEP::Hep3Vector     CofM;
      int                   n_seeds;
      int                   fFirstStation;
      int                   fLastStation;
      McPart_t*             fMcPart;
      int                   fNHits;
      int                   fNHitsMcP;	// Nhits by the "best" particle"
      LsqSums2              fTzSums;

      DeltaCandidate() {
	fNumber = -1;
	for(int s=0; s<kNStations; ++s) {
	  st_used[s] = false;
	  dxy    [s] = -1;
	  seed   [s] = NULL;
	}
	fFirstStation = -1;
	fLastStation  = -1;
	fMcPart       = NULL;
	fNHits        = 0;
	fTzSums.clear();
      }

      double PredictedTime(double Z) { return fTzSums.yMean()+fTzSums.dydx()*(Z-CofM.z()); }
    };
//-----------------------------------------------------------------------------
// data structure passed to the histogramming routine
//-----------------------------------------------------------------------------
    struct Data_t {
      const art::Event*            event;
      const TTracker*              tracker;
      std::string                  strawDigiMCCollectionTag;
      std::string                  ptrStepPointMCVectorCollectionTag;
      std::vector<DeltaSeed*>      seedHolder [kNStations];
      std::vector<DeltaCandidate>  deltaCandidateHolder;
      PanelZ_t                     oTracker[kNStations][kNFaces][kNPanelsPerFace];
      int                          nseeds;
      int                          nseeds_per_station[kNStations];
      const StrawHitCollection*    shcol;
      int                          nsh; // total number of straw hits
      int                          testOrder;         // 
      int                          printElectrons;    //
      int                          printElectronsMinNHits;
      float                        printElectronsMaxFReco;
      float                        printElectronsMinMom;
      int                          printDeltaSeeds;
      int                          printDeltaCandidates;
      int                          debugLevel;	     // printout level
      int                          diagLevel;            // histogramming level
      bool                         mcDiag;
      int                          printOTracker;     //
    };
  }
}
#endif
