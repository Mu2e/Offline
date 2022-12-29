#ifndef __CalPatRec_DeltaFinder_types_hh__
#define __CalPatRec_DeltaFinder_types_hh__

namespace art {
  class Event;
};

namespace fhicl {
  class ParameterSet;
};

#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/RecoDataProducts/inc/StereoHit.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

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
      double     x;                        // x-coordinate of the intersection point
      double     y;                        // y-coordinate of the intersection point
      double     z;                        // y-coordinate of the intersection point
      double     wd1;                      // distance btw the 1st hit and the intersection
      double     wd2;                      // distance btw the 2nd hit and the intersection
    };

    struct DeltaCandidate;

    enum {
      kNStations      = StrawId::_nplanes/2,   // number of tracking stations
      kNFaces         = StrawId::_nfaces*2 ,   // N(faces) per station (4)
      kNPanelsPerFace = StrawId::_npanels/2    // = 3
    };

    struct HitData_t {
      const ComboHit*         fHit;
      int                     fSeedIndex;
      int                     fNSecondHits;
      int                     fDeltaIndex;
      float                   fChi2Min;
      float                   fSigW;     // cached resolution along the wire
      float                   fDr;        // work variable

      HitData_t(const ComboHit* Hit, float SigW) {
        fHit         = Hit;
        fChi2Min     = 999999.;
        fSigW        = SigW;
        fSeedIndex   = -1;
        fNSecondHits =  0;
        fDeltaIndex  = -1;
        fDr          = 1.1e10;
      }

      int Used() const { return (fSeedIndex >= 0) ; }
    };

    struct PanelZ_t {
      int                              fNHits  ; // guess, total number of ComboHits
      std::vector<HitData_t>           fHitData;
      const Panel*                     fPanel;      // backward pointer to the tracker panel
      double                           wx;          // direction cosines of the wires, assumed to be all the same
      double                           wy;
      double                           phi;         // phi angle of the wire
      double                           z;           //
      float                            tmin;        // for hits stored on this panel
      float                            tmax;
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
      int                           fID;         // SimParticle::id().asInt()
      int                           fPdgID;
      int                           fNHitsCE;
      int                           fNHitsDelta; // number of hits associated with reconstructed delta electrons
      float                         fTime;
      float                         fHitTMin;    // min and max times of the straw hits
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
        fTime         = 1.e6;                 // at initialization, make it absurd
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
      int                            fNumber;           // number within the station
      int                            fStation;          // station with seed stereo hit
      int                            fType;             // defines indices of the two faces used for preseed seach
      int                            fGood;             // <-killer number> if not to be used - what about 0 ?
      int                            fNFacesWithHits;
      int                            fNHits;
      int                            fNStrawHits;       // total number of hits
      int                            fNHitsCE;

      int                            fSFace[2];         // faces making the stereo seed
      // HitData_t*                     fHitData[2];    // stereo hit seeding the seed search
      float                          fChi21;            // chi2's of the two initial hits, also stored in hit data
      float                          fChi22;
                                                        // 0: used in recovery
      int                            fFaceProcessed[kNFaces];
      HitData_t*                     hitData       [kNFaces];
      McPart_t*                      fMcPart       [kNFaces];
      XYZVectorF                     CofM;
      float                          fPhi;              // cache to speed up the phi checks
      float                          fMinHitTime;          // min and max times of the included hits
      float                          fMaxHitTime;
      McPart_t*                      fPreSeedMcPart[2]; // McPart_t's corresponding to initial intersection
      int                            fDeltaIndex;
      float                          fChi2All;

      DeltaSeed ();
      DeltaSeed (int Index, int Station, int Face0, HitData_t* Hd0, int Face1, HitData_t* Hd1);
      ~DeltaSeed() {}

      int              SFace(int I)       { return fSFace[I]; }
      float            Chi2N   ()         { return (fChi21+fChi22)/fNHits; }
      float            Chi2Tot ()         { return (fChi21+fChi22); }
      float            Chi2All ()         { return fChi2All; }
      float            Chi2AllDof ()      { return fChi2All/fNHits; }
      HitData_t*       HitData (int Face) { return hitData[Face]; } // no boundary check !
      int              NHits   ()         { return fNHits; }
      int              NStrawHits()       { return fNStrawHits; }
      int              MCTruth ()         { return (fPreSeedMcPart[0] != NULL) && (fPreSeedMcPart[0] == fPreSeedMcPart[1]) ; }
      bool             Used    ()         { return (fDeltaIndex >= 0); }
      int              Good    ()         { return (fGood       >= 0); }

      float            Phi     () const   { return fPhi; }
//-----------------------------------------------------------------------------
// drift time can't be < 0
// fMaxTime < particle T0 < fMinTime
//-----------------------------------------------------------------------------
      float            Time      () { return (fMinHitTime+fMaxHitTime)/2;    }
      float            MinHitTime() { return fMinHitTime; }
      float            MaxHitTime() { return fMaxHitTime; }
//-----------------------------------------------------------------------------
// assumed range of allowed hit times for this seed
//-----------------------------------------------------------------------------
      float            T0Min     () { return (fMinHitTime+fMaxHitTime)/2-25; }
      float            T0Max     () { return (fMinHitTime+fMaxHitTime)/2+25; }
    };

    struct DeltaCandidate {
    public:
      int                   fIndex;                 // >= 0: index, <0: -1000-index merged
      int                   fFirstStation;
      int                   fLastStation;
      DeltaSeed*            seed   [kNStations];
      float                 dxy    [kNStations];   // used only for diagnostics
      float                 fT0Min [kNStations];   // acceptable hit times (no need to account for the drift time!)
      float                 fT0Max [kNStations];
      XYZVectorF            CofM;
      float                 phi;
      int                   n_seeds;
      McPart_t*             fMcPart;
      int                   fNHits;
      int                   fNStrawHits;
      int                   fNHitsMcP;               // Nhits by the "best" particle"
      int                   fNHitsCE;

      DeltaCandidate();
      DeltaCandidate(int Index, DeltaSeed* Seed, int Station);

      int        Active               () const { return (fIndex >= 0); }
      int        Index                () const { return fIndex ; }
      int        NSeeds               () { return n_seeds; }
      int        NHits                () { return fNHits; }
      int        NStrawHits           () { return fNStrawHits; }
      DeltaSeed* Seed            (int I) { return seed[I]; }
      bool       StationUsed     (int I) { return (seed[I] != NULL); }
      float      T0Min           (int I) { return fT0Min[I]; }
      float      T0Max           (int I) { return fT0Max[I]; }
      float      Time            (int I) { return (fT0Max[I]+fT0Min[I])/2.; }
      int        LastStation          () { return fLastStation ; }
      int        FirstStation         () { return fFirstStation; }

      void       AddSeed            (DeltaSeed* Ds, int Station);

      void       MergeDeltaCandidate(DeltaCandidate* Delta);

      void       SetIndex(int Index) { fIndex = Index; }
    };
//-----------------------------------------------------------------------------
// data structure passed to the diagnostics plugin
//-----------------------------------------------------------------------------
    struct Data_t {
      const art::Event*             event;
      const Tracker*                tracker;
      art::InputTag                 chCollTag;
      art::InputTag                 chfCollTag;

      art::InputTag                 sdmcCollTag;
      std::vector<DeltaSeed*>       listOfSeeds[kNStations]; // seeds with the first station being this
      std::vector<DeltaCandidate>   listOfDeltaCandidates;

      PanelZ_t                      oTracker   [kNStations][kNFaces][kNPanelsPerFace];
      int                           stationUsed[kNStations];

      int                           nseeds;
      int                           nseeds_per_station[kNStations];
      const ComboHitCollection*     chcol;
      StrawHitFlagCollection*       chfcol;       // output combohit flags
      const TimeClusterCollection*  tpeakcol;
      int                           debugLevel;   // printout level

      DeltaCandidate*     deltaCandidate(int I) { return &listOfDeltaCandidates[I]; }
      DeltaSeed*          deltaSeed(int Station, int I) { return listOfSeeds[Station][I]; }
    };

//-----------------------------------------------------------------------------
// finally, utility functions
//-----------------------------------------------------------------------------
    int findIntersection(const HitData_t* Hd1, const HitData_t* Hd2, Intersection_t* Result);
  }
}
#endif
