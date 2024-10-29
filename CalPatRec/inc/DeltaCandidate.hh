#ifndef CalPatRec_DeltaCandidate_hh
#define CalPatRec_DeltaCandidate_hh

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
      DeltaSeed*            fSeed   [kNStations];
      XYZVectorF            CofM;
      float                 fNx;                   //
      float                 fNy;                   //
      int                   fNSeeds;
      int                   fNHits;                // n(combo hits)
      int                   fNStrawHits;           // number of straw hits
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
      DeltaCandidate(int Index, DeltaSeed* Seed);

      int        Active               () const { return (fIndex >= 0); }
      int        Index                () const { return fIndex ; }
      int        Mask                 () const { return fMask  ; }
      int        NSeeds               () const { return fNSeeds; }
      int        nHits                () const { return fNHits; }
      int        nStrawHits           () const { return fNStrawHits; }
      DeltaSeed* Seed            (int I) const { return fSeed[I]; }
      bool       StationUsed     (int I) const { return (fSeed[I] != NULL); }
      int        LastStation          () const { return fLastStation ; }
      int        FirstStation         () const { return fFirstStation; }
      float      EDep                 () const { return fSumEDep/fNStrawHits; }
      double     Xc                   () const { return CofM.x(); }
      double     Yc                   () const { return CofM.y(); }
      double     Rho                  () const { return CofM.Rho(); }
      double     Nx                   () const { return fNx ; }
      double     Ny                   () const { return fNy ; }
//-----------------------------------------------------------------------------
// add Seed  (station is defined by the Seed
//-----------------------------------------------------------------------------
      void       AddSeed            (DeltaSeed*      Ds);
//-----------------------------------------------------------------------------
// diagnostics: fill 'this' with the parameters of 'Delta' w/o segment in a
//              given 'Station'
//-----------------------------------------------------------------------------
      void       removeSeed         (const DeltaCandidate* Delta, int Station);
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
