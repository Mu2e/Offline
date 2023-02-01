#ifndef __CalPatRec_DeltaSeed_hh__
#define __CalPatRec_DeltaSeed_hh__

#include "TObject.h"
#include "TClonesArray.h"
#include "Math/Vector3D.h"
// #include "Math/Vector4D.h"

using ROOT::Math::XYZVectorF;

#include "Offline/CalPatRec/inc/DeltaFinder_enums.hh"

namespace mu2e {

  namespace DeltaFinderTypes {
    struct HitData_t;
    struct McPart_t;
  }

  using  DeltaFinderTypes::HitData_t;
  using  DeltaFinderTypes::McPart_t;

  class  Panel;
  class  SimParticle;

  class DeltaSeed : public TObject {
  public:
    int          fIndex;            // index within the station
    int          fStation;          // station with seed stereo hit
    int          fType;             // defines indices of the two faces used for preseed seach
    int          fGood;             // <killer number> if not to be used - what about 0 ?
    int          fNFacesWithHits;
    int          fNHits;            // total number of combo hits
    int          fNStrawHits;       // total number of straw hits
    int          fNHitsCE;          // number of associated CE hits

    int          fSFace[2];         // faces making the stereo seed
    float        fChi21;            // chi2's of the two initial hits, also stored in hit data
    float        fChi22;
                                      // 0: used in recovery
    int          fFaceProcessed[kNFaces];

    HitData_t*   hitData       [kNFaces];

    McPart_t*    fMcPart       [kNFaces];
                                      // XY coordinate sums
    double       fSnx2;
    double       fSnxy;
    double       fSny2;
    double       fSnxr;
    double       fSnyr;
    float        fSumEDep;             // sum over the straw hits

    XYZVectorF   CofM;                 // COG
    //    float        fPhi;                 // cache to speed up the phi checks
    float        fZ;                   // Z-coordinate of the center of the corresponding station

    double       fSumT;

    float        fMinHitTime;          // min and max times of the included hits
    float        fMaxHitTime;
    McPart_t*    fPreSeedMcPart[2];    // McPart_t's for initial intersection
    int          fDeltaIndex;
                                       // chi2's
    float        fChi2All;
    float        fChi2Perp;
    float        fChi2Delta;           // chi2 when the seed is added to Delta

    DeltaSeed () {}
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
    float            Chi2Delta()        { return fChi2Delta   ; }
    HitData_t*       HitData (int Face) { return hitData[Face]; } // no boundary check !
    int              NHits   ()         { return fNHits; }
    int              NHitsCE ()         { return fNHitsCE; }
    int              NStrawHits()       { return fNStrawHits; }
    float            SumEDep ()         { return fSumEDep ; }
    float            EDep    ()         { return fSumEDep/fNStrawHits ; }
    int              MCTruth ()         { return (fPreSeedMcPart[0] != NULL) && (fPreSeedMcPart[0] == fPreSeedMcPart[1]) ; }
    bool             Used    ()         { return (fDeltaIndex >= 0); }
    int              Good    ()         { return (fGood       >= 0); }
    float            TMean   ()         { return fSumT/fNStrawHits;  }
    float            Z       ()         { return fZ;  }
    double           Xc      () const { return CofM.x(); }
    double           Yc      () const { return CofM.y(); }
    double           R       () const { return CofM.R(); }

    //    float            Phi     () const   { return fPhi; }
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
    void             Chi2(double Xc, double Yc, double SigmaR2, double& Chi2All, double& Chi2Perp);

    ClassDef(mu2e::DeltaSeed,0)
  };
}
#endif
