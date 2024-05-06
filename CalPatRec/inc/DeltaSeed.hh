#ifndef CalPatRec_DeltaSeed_hh
#define CalPatRec_DeltaSeed_hh

#include "TObject.h"
#include "TClonesArray.h"
#include "Math/Vector3D.h"

using ROOT::Math::XYZVectorF;

#include "Offline/CalPatRec/inc/DeltaFinder_structures.hh"

namespace mu2e {

  namespace DeltaFinderTypes {
    struct McPart_t;
  }

  using  CalPatRec::HitData_t;
  using  DeltaFinderTypes::McPart_t;

  class  Panel;
  class  SimParticle;

  class DeltaSeed {
  public:
    static       float fSigT2;      // resolution in Tcorr, per hit (8 ns^2)
  public:
    int          fIndex;            // seed index within the station
    int          fStation;          // delta seed station
    int          fType;             // defines indices of the two faces used for preseed seach
    int          fGood;             // <killer number> if not to be used - what about 0 ?
    int          fNHits;            // total number of combo hits
    int          fNStrawHits;       // total number of straw hits

    int          fSFace[2];         // faces making the stereo seed
    float        fChi21;            // chi2's of the two initial hits, also stored in hit data
    float        fChi22;
                                      // 0: used in recovery
    int          fFaceProcessed[kNFaces];

    HitData_t*   fHitData      [kNFaces];

                                        // XY coordinate sums
    double       fSnx2;
    double       fSnxy;
    double       fSny2;
    double       fSnxr;
    double       fSnyr;
    float        fSumEDep;             // sum over the straw hits

    XYZVectorF   CofM;                 // COG
    float        fZ;                   // Z-coordinate of the center of the corresponding station

    float        fSumT;
    float        fSumT2;

    float        fMinHitTime;          // min and max times of the included hits
    float        fMaxHitTime;
    int          fDeltaIndex;
    int          fProtonIndex;
                                       // chi2's
    float        fChi2Par;
    float        fChi2Perp;
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
    DeltaSeed (int Index) { fIndex = Index; }

    ~DeltaSeed() {}

    void             Init(HitData_t* Hd0, HitData_t* Hd1, float Xc, float Yx, float Zc);

    int              Station ()         { return fStation; }
    int              Index   ()         { return fIndex; }
    int              SFace   (int I)    { return fSFace[I]; }

    float            Chi2Par  ()        { return fChi2Par; }
    float            Chi2Perp ()        { return fChi2Perp; }
    float            Chi2Tot  ()        { return (fChi2Par+fChi2Perp); }
    float            Chi2ParN ()        { return fChi2Par /fNHits; }
    float            Chi2PerpN()        { return fChi2Perp/fNHits; }
    float            Chi2TotN ()        { return (fChi2Par+fChi2Perp)/fNHits; }
//-----------------------------------------------------------------------------
// 10000+idelta means that all hits of this segments, one by one, have been
// picked up by the same delta candidate 'idelta'
//-----------------------------------------------------------------------------
    int              deltaIndex ()      { return fDeltaIndex  % 10000 ; }
    int              protonIndex()      { return fProtonIndex ; }

    HitData_t*       HitData (int Face) { return fHitData[Face]; } // no boundary check !
    int              nHits   ()         { return fNHits; }
    int              nStrawHits()       { return fNStrawHits; }

    float            SumEDep ()         { return fSumEDep ; }
    float            EDep    ()         { return fSumEDep/fNStrawHits ; }
    bool             Used    ()         { return (fDeltaIndex >= 0); }
    int              Good    ()         { return (fGood       >= 0); }

    float            TMean   ()         { return fSumT/fNHits;  }
    float            Chi2Time ()        { return (fSumT2/fNHits-(fSumT/fNHits)*(fSumT/fNHits))/fSigT2; }

    float            Z       ()         { return fZ;  }
    double           Xc      () const   { return CofM.x(); }
    double           Yc      () const   { return CofM.y(); }
    double           Rho     () const   { return CofM.Rho(); }
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

    float            Dt12      () {
      return fHitData[fSFace[0]]->fHit->time()-fHitData[fSFace[1]]->fHit->time();
    }

    float            DtCorr12      () {
      return fHitData[fSFace[0]]->fCorrTime-fHitData[fSFace[1]]->fCorrTime;
    }

    void             SetStation    (int Station) { fStation     = Station; }
    void             SetDeltaIndex (int Index  ) { fDeltaIndex  = Index  ; }
    void             setProtonIndex(int Index  ) { fProtonIndex = Index  ; }
//-----------------------------------------------------------------------------
// less trivial functions .. HitData_t knows its ZFace
//-----------------------------------------------------------------------------
    void             AddHit             (HitData_t* Hd);
    void             ReplaceFirstHit    (HitData_t* Hd);
    void             CalculateCogAndChi2(float RCore, float SigmaR2);
    void             Chi2(float Xc, float Yc, float RCore, float SigmaR2, float& Chi2All, float& Chi2Perp);
  };
}
#endif
