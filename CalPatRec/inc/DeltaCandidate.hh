#ifndef __CalPatRec_DeltaCandidate_hh__
#define __CalPatRec_DeltaCandidate_hh__

#include "Offline/CalPatRec/inc/DeltaFinder_structures.hh"
#include "Offline/CalPatRec/inc/DeltaSeed.hh"

namespace mu2e {
  class Panel;
  class SimParticle;

    struct DeltaCandidate {
      enum  {
        kEDepBit   = 0x00000001, // <<  0
        kNHitsBit  = 0x00000002
      };

    public:
      int                   fIndex;                 // >= 0: index, <0: -1000-index merged
      int                   fMask;                  // bitmask , if zero, the candiate is good
      int                   fFirstStation;
      int                   fLastStation;
      DeltaSeed*            seed   [kNStations];
      float                 dxy    [kNStations];   // used only for diagnostics
      // float                 fT0Min [kNStations];   // acceptable hit times (no need to account for the drift time!)
      // float                 fT0Max [kNStations];
      XYZVectorF            CofM;
      int                   fNSeeds;
      McPart_t*             fMcPart;               // "best" MC particle
      int                   fNHits;                // n(combo hits)
      int                   fNStrawHits;           // number of straw hits
      int                   fNHitsMcP;             // N combo hits by the "best" particle"
      int                   fNHitsCE;              // N(hits) by CE
      float                 fSumEDep;              //
                                                   // LSQ sums
      double                fSx;
      double                fSy;
      double                fSnx2;
      double                fSnxy;
      double                fSny2;
      double                fSnxr;
      double                fSnyr;
                                                   // LSQ sumz for TZ
      double                fSt;
      double                fSz;
      double                fSt2;
      double                fStz;
      double                fSz2;

      double                fT0;
      double                fDtDz;
      double                fSigT0;
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
      DeltaCandidate();
      DeltaCandidate(int Index, DeltaSeed* Seed, int Station);

      int        Active               () const { return (fIndex >= 0); }
      int        Index                () const { return fIndex ; }
      int        Mask                 () const { return fMask  ; }
      int        NSeeds               () const { return fNSeeds; }
      int        NHits                () const { return fNHits; }
      int        NHitsMcP             () const { return fNHitsMcP; }
      int        NStrawHits           () const { return fNStrawHits; }
      DeltaSeed* Seed            (int I) const { return seed[I]; }
      bool       StationUsed     (int I) const { return (seed[I] != NULL); }
      // float      T0Min           (int I) const { return fT0Min[I]; }
      // float      T0Max           (int I) const { return fT0Max[I]; }
      // float      Time            (int I) const { return (fT0Max[I]+fT0Min[I])/2.; }
      int        LastStation          () const { return fLastStation ; }
      int        FirstStation         () const { return fFirstStation; }
      float      EDep                 () const { return fSumEDep/fNStrawHits; }
      float      FBest                () const { return float(fNHitsMcP)/fNHits; }
      double     Xc                   () const { return CofM.x(); }
      double     Yc                   () const { return CofM.y(); }
      double     R                    () const { return CofM.R(); }
      double     Nx                   () const { return CofM.x()/CofM.R() ; }
      double     Ny                   () const { return CofM.y()/CofM.R() ; }

      void       AddSeed            (DeltaSeed*      Ds   , int Station);
      void       markHitsAsUsed     ();
      void       MergeDeltaCandidate(DeltaCandidate* Delta, int PrintErrors);

      void       SetIndex(int Index) { fIndex = Index; }

      float      T0(int Station)     {
        if (fNSeeds == 1) return fT0;
        else              return fT0 + fDtDz*DeltaFinderTypes::stationZ[Station];
      }
    };
}
#endif
