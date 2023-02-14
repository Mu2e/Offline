#ifndef __CalPatRec_ProtonCandidate_hh__
#define __CalPatRec_ProtonCandidate_hh__

#include "Offline/CalPatRec/inc/DeltaFinder_structures.hh"
#include "Offline/CalPatRec/inc/DeltaSeed.hh"

namespace mu2e {
  class Panel;
  class SimParticle;

    struct ProtonCandidate {
      enum  {
        kEDepBit   = 0x00000001, // <<  0
        kNHitsBit  = 0x00000002
      };

    public:
      int                     fIndex;                 // >= 0: index, <0: -1000-index merged
      int                     fMask;                  // bitmask , if zero, the candiate is good
      int                     fFirstStation;
      int                     fLastStation;
      std::vector<HitData_t*> fHitData[kNStations][kNFaces];
      int                     fPanelID[kNStations][kNFaces];
      int                     fNHitsStation[kNStations];
      float                   fMean        [kNStations];

      int                     fNStationsWithHits;
      int                     fNHits;              // total N(combo hits)
      int                     fNStrawHits;         // total N(straw hits)

                                                   // LSQ sumz for TZ
      double                  fSt;
      double                  fSz;
      double                  fSt2;
      double                  fStz;
      double                  fSz2;

      float                   fSumEDep;

      double                  fT0;
      double                  fDtDz;
      double                  fSigT0;

      McPart_t*               fMcPart;             // "best" MC particle
      int                     fNHitsMcP;           // N combo hits by the "best" particle"
      int                     fNHitsCE;            // N(hits) by CE
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
      ProtonCandidate(int Index);

      int        Active               () const { return (fIndex >= 0); }
      int        Index                () const { return fIndex ; }
      int        Mask                 () const { return fMask  ; }
      int        NStationsWithHits    () const { return fNStationsWithHits; }
      int        NHitsTotal           () const { return fNHits; }

      int        nHits     (int Station, int Face) const  { return fHitData[Station][Face].size(); }
      HitData_t* hitData   (int Station, int Face, int I) { return fHitData[Station][Face][I]; }


      int        NHitsMcP             () const { return fNHitsMcP; }
      int        NStrawHits           () const { return fNStrawHits; }
      // float      T0Min           (int I) const { return fT0Min[I]; }
      // float      T0Max           (int I) const { return fT0Max[I]; }
      int        LastStation          () const { return fLastStation ; }
      int        FirstStation         () const { return fFirstStation; }
      float      EDep                 () const { return fSumEDep/fNStrawHits; }
      float      FBest                () const { return float(fNHitsMcP)/fNHits; }
//-----------------------------------------------------------------------------
// add Seed  (station is defined by the Seed
//-----------------------------------------------------------------------------
      void       addSeed            (DeltaSeed* Seed);
      void       addHit             (int Station, HitData_t* Hit, int UpdateTime = 1);
//-----------------------------------------------------------------------------
// in case of one station, fDtDz = 0
//-----------------------------------------------------------------------------
      float      T0(int Station)     { return fT0 + fDtDz*DeltaFinderTypes::stationZ[Station]; }

//-----------------------------------------------------------------------------
// an attempt to predict phi in a given station, returns -100 if no prediction
//-----------------------------------------------------------------------------
      float      Phi(int Station);
    };
}
#endif
