#ifndef __CalPatRec_DeltaFinder_types_hh__
#define __CalPatRec_DeltaFinder_types_hh__

namespace art {
  class Event;
};

namespace fhicl {
  class ParameterSet;
};

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/RecoDataProducts/inc/StereoHit.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

#include "Offline/Mu2eUtilities/inc/McUtilsToolBase.hh"

namespace mu2e {
  class Panel;
  class SimParticle;
//-----------------------------------------------------------------------------
// delta-electron seed: structure within the station
// doesn't own anything, no need to delete any pinters
//-----------------------------------------------------------------------------
  namespace DeltaFinderTypes {

    struct Config {
      fhicl::Atom<std::string> tool_type             {fhicl::Name("tool_type"             ), fhicl::Comment("tool type: DeltaFinderDiag")     };
      // fhicl::Atom<art::InputTag> spmcCollTag         {fhicl::Name("spmcCollTag"           ), fhicl::Comment("StepPointMC coll tag"      )     };
      fhicl::Atom<int>         mcTruth               {fhicl::Name("mcTruth"               ), fhicl::Comment("MC truth")                       };
      fhicl::Atom<int>         diagLevel             {fhicl::Name("diagLevel"             ), fhicl::Comment("diagnostic level")               };
      fhicl::Atom<bool>        mcDiag                {fhicl::Name("mcDiag"                ), fhicl::Comment("MC diag")                        };
      fhicl::Atom<int>         printOTracker         {fhicl::Name("printOTracker"         ), fhicl::Comment("print ordered Tracker")          };
      fhicl::Atom<int>         printComboHits        {fhicl::Name("printComboHits"        ), fhicl::Comment("print combo hits")               };
      fhicl::Atom<int>         printElectrons        {fhicl::Name("printElectrons"        ), fhicl::Comment("print electrons")                };
      fhicl::Atom<int>         printElectronsHits    {fhicl::Name("printElectronsHits"    ), fhicl::Comment("print electron hits")            };
      fhicl::Atom<int>         printElectronsMinNHits{fhicl::Name("printElectronsMinNHits"), fhicl::Comment("minNhhits for printed electrons")};
      fhicl::Atom<float>       printElectronsMaxFReco{fhicl::Name("printElectronsMaxFReco"), fhicl::Comment("maxFReco for printed electrons" )};
      fhicl::Atom<float>       printElectronsMinMom  {fhicl::Name("printElectronsMinMom"  ), fhicl::Comment("min mom for printed electrons"  )};
      fhicl::Atom<float>       printElectronsMaxMom  {fhicl::Name("printElectronsMaxMom"  ), fhicl::Comment("max mom for printed electrons"  )};
      fhicl::Atom<int>         printDeltaSeeds       {fhicl::Name("printDeltaSeeds"       ), fhicl::Comment("if 1, print delta seeds"        )};
      fhicl::Atom<int>         printDeltaCandidates  {fhicl::Name("printDeltaCandidates"  ), fhicl::Comment("if 1, print delta candidates"   )};
      fhicl::Atom<int>         printShcol            {fhicl::Name("printShcol"            ), fhicl::Comment("if 1, print shColl"             )};

      fhicl::Table<McUtilsToolBase::Config> mcUtils{fhicl::Name("mcUtils"       ), fhicl::Comment("MC Diag plugin") };
    };
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
      float                   fSigW2;    // cached resolution^2 along the wire
      float                   fCorrTime; // cache hit corrected time

      HitData_t(const ComboHit* Hit) {
        fHit         = Hit;
        fChi2Min     = 99999.0;
        float sigw   =  Hit->posRes(ComboHit::wire);
        fSigW2       = sigw*sigw;  // to be used to calculate chi2...
        fSeedIndex   = -1;
        fNSecondHits =  0;
        fDeltaIndex  = -1;
        fCorrTime    = Hit->correctedTime();
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
      int                            fIndex;            // index within the station
      int                            fStation;          // station with seed stereo hit
      int                            fType;             // defines indices of the two faces used for preseed seach
      int                            fGood;             // <killer number> if not to be used - what about 0 ?
      int                            fNFacesWithHits;
      int                            fNHits;            // total number of combo hits
      int                            fNStrawHits;       // total number of straw hits
      int                            fNHitsCE;          // number of associated CE hits

      int                            fSFace[2];         // faces making the stereo seed
      float                          fChi21;            // chi2's of the two initial hits, also stored in hit data
      float                          fChi22;
                                                        // 0: used in recovery
      int                            fFaceProcessed[kNFaces];
      HitData_t*                     hitData       [kNFaces];
      McPart_t*                      fMcPart       [kNFaces];
                                                        // XY coordinate sums
      double                         fSx;
      double                         fSy;
      double                         fSnx2;
      double                         fSnxny;
      double                         fSny2;
      double                         fSnxnr;
      double                         fSnynr;
      float                          fSumEDep;          // sum over the straw hits

      XYZVectorF                     CofM;                 // COG
      float                          fPhi;                 // cache to speed up the phi checks
      float                          fMinHitTime;          // min and max times of the included hits
      float                          fMaxHitTime;
      McPart_t*                      fPreSeedMcPart[2];    // McPart_t's for initial intersection
      int                            fDeltaIndex;
                                                           // chi2's
      float                          fChi2All;
      float                          fChi2Perp;

      DeltaSeed ();
      DeltaSeed (int Index, int Station, int Face0, HitData_t* Hd0, int Face1, HitData_t* Hd1);
      ~DeltaSeed() {}

      int              Station ()         { return fStation; }
      int              Index   ()         { return fIndex; }
      int              SFace(int I)       { return fSFace[I]; }
      float            Chi2TotN()         { return (fChi21+fChi22)/fNHits; }
      float            Chi2Tot ()         { return (fChi21+fChi22); }
      float            Chi2All ()         { return fChi2All; }
      float            Chi2Perp()         { return fChi2Perp; }
      float            Chi2PerpN()        { return fChi2Perp/fNHits; }
      float            Chi2AllN ()        { return fChi2All/fNHits ; }
      HitData_t*       HitData (int Face) { return hitData[Face]; } // no boundary check !
      int              NHits   ()         { return fNHits; }
      int              NHitsCE ()         { return fNHitsCE; }
      int              NStrawHits()       { return fNStrawHits; }
      float            SumEDep ()         { return fSumEDep ; }
      float            EDep    ()         { return fSumEDep/fNStrawHits ; }
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
      float            T0Min     () { return (fMinHitTime+fMaxHitTime)/2-20; }
      float            T0Max     () { return (fMinHitTime+fMaxHitTime)/2+20; }
//-----------------------------------------------------------------------------
// less trivial functions
//-----------------------------------------------------------------------------
      void             AddHit         (HitData_t* Hd, int Face);
      void             ReplaceFirstHit(HitData_t* Hd);
      void             CalculateCogAndChi2(double SigmaR2);
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
      int                   fNHits;     // n(combo hits)
      int                   fNStrawHits;
      int                   fNHitsMcP;               // N combo hits by the "best" particle"
      int                   fNHitsCE;
      float                 fSumEDep;      //

      DeltaCandidate();
      DeltaCandidate(int Index, DeltaSeed* Seed, int Station);

      int        Active               () const { return (fIndex >= 0); }
      int        Index                () const { return fIndex ; }
      int        NSeeds               () { return n_seeds; }
      int        NHits                () { return fNHits; }
      int        NHitsMcP             () { return fNHitsMcP; }
      int        NStrawHits           () { return fNStrawHits; }
      DeltaSeed* Seed            (int I) { return seed[I]; }
      bool       StationUsed     (int I) { return (seed[I] != NULL); }
      float      T0Min           (int I) { return fT0Min[I]; }
      float      T0Max           (int I) { return fT0Max[I]; }
      float      Time            (int I) { return (fT0Max[I]+fT0Min[I])/2.; }
      int        LastStation          () { return fLastStation ; }
      int        FirstStation         () { return fFirstStation; }
      float      EDep                 () { return fSumEDep/fNStrawHits; }
      float      FBest                () { return float(fNHitsMcP)/fNHits; }

      void       AddSeed            (DeltaSeed*      Ds   , int Station);
      void       MergeDeltaCandidate(DeltaCandidate* Delta, int PrintErrors);

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

      DeltaCandidate* deltaCandidate(int I)              { return &listOfDeltaCandidates[I]; }
      DeltaSeed*      deltaSeed     (int Station, int I) { return listOfSeeds[Station][I]; }

      void printHitData       (HitData_t*      HitData, const char* Option = "");
      void printDeltaSeed     (DeltaSeed*      Seed   , const char* Option = "");
      void printDeltaCandidate(DeltaCandidate* Delta  , const char* Option = "");
    };

//-----------------------------------------------------------------------------
// finally, utility functions
//-----------------------------------------------------------------------------
    int findIntersection(const HitData_t* Hd1, const HitData_t* Hd2, Intersection_t* Result);
  }
}
#endif
