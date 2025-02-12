#ifndef CalPatRec_ProtonCandidate_hh
#define CalPatRec_ProtonCandidate_hh

#include "Offline/CalPatRec/inc/DeltaFinder_structures.hh"
#include "Offline/CalPatRec/inc/DeltaSeed.hh"

namespace mu2e {
  class Panel;
  class SimParticle;

  using DeltaFinderTypes::PhiPrediction_t;

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
      std::vector<HitData_t*> fHitData     [kNStations][kNFaces];
      int                     fPanelID     [kNStations][kNFaces];
      int                     fNHitsStation[kNStations];
      int                     fMinHitTime  [kNStations];
      int                     fMaxHitTime  [kNStations];
      float                   fSumX        [kNStations];
      float                   fSumY        [kNStations];
      float                   fPhi         [kNStations];

      int                     fNStationsWithHits;
      int                     fNHitsTot;           // total N(combo hits)
      int                     fNStrawHitsTot;      // total N(straw hits)

                                                   // LSQ sumz for TZ
      double                  fSt;
      double                  fSz;
      double                  fSt2;
      double                  fStz;
      double                  fSz2;

      float                   fTMid;               // time in the station closest to the center

      float                   fSumEDep;

      float                   fT0;
      float                   fDtDz;
      float                   fSigT0;

      McPart_t*               fMcPart;             // "best" MC particle
      int                     fNHitsMcP;           // N combo hits by the "best" particle"
      int                     fNHitsCE;            // N(hits) by CE
      int                     fTimeIndex;
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
      ProtonCandidate(int Index);
      ProtonCandidate() : ProtonCandidate(0) {}

      void       init                 ();

      int        Active               () const { return (fIndex >= 0); }
      int        index                () const { return fIndex ; }
      int        Mask                 () const { return fMask  ; }
      int        nStationsWithHits    () const { return fNStationsWithHits; }
      int        nHitsTot             () const { return fNHitsTot; }

      int        nHits     (int Station, int Face) const  { return fHitData[Station][Face].size(); }
      HitData_t* hitData   (int Station, int Face, int I) { return fHitData[Station][Face][I]; }

      int        nHitsStation(int Station) const  { return fNHitsStation[Station]; }


      int        NHitsMcP             () const { return fNHitsMcP; }
      int        nStrawHitsTot        () const { return fNStrawHitsTot; }

      float      minHitTime(int Station) const { return fMinHitTime[Station]; }
      float      maxHitTime(int Station) const { return fMaxHitTime[Station]; }
      int        LastStation          () const { return fLastStation ; }
      int        FirstStation         () const { return fFirstStation; }
      float      eDep                 () const { return fSumEDep/fNStrawHitsTot; }
      float      FBest                () const { return float(fNHitsMcP)/fNHitsTot; }

      float      xMean(int Station)  { return fSumX[Station]/fNHitsStation[Station]; }
      float      yMean(int Station)  { return fSumY[Station]/fNHitsStation[Station]; }
      float      phi  (int Station)  { return fPhi[Station] ; }
      float      t0   (int Station)  { return fT0 + fDtDz*DeltaFinderTypes::stationZ[Station]; }
      float      tMid()              { return fTMid; }
      int        timeIndex()         { return fTimeIndex; }
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
      void       setIndex(int Index) { fIndex = Index; }
//-----------------------------------------------------------------------------
// add Seed  (station is defined by the Seed
//-----------------------------------------------------------------------------
      void       addSeed            (DeltaSeed* Seed);
      void       addHit             (int Station, HitData_t* Hit, int UpdateTime = 1);
      void       removeHit          (int Station, HitData_t* Hit, int UpdateTime = 1);
      void       updateTime         ();
//-----------------------------------------------------------------------------
// in case of one station, fDtDz = 0
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// an attempt to predict phi in a given station, returns -100 if no prediction
// if the panel ID can be predicted, use that
//-----------------------------------------------------------------------------
      void       predictPhi(int Station, PhiPrediction_t* Prediction);
    };
}
#endif
