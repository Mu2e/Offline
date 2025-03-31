///////////////////////////////////////////////////////////////////////////////
// ComboHit.hh needs the definition of ProductID ...
///////////////////////////////////////////////////////////////////////////////
#ifndef CalPatRec_MCPart_t_hh
#define CalPatRec_MCPart_t_hh

#include "TObject.h"

namespace mu2e {

  class  Panel;
  class  SimParticle;
  struct DeltaCandidate;
  struct ProtonCandidate;
  class  DeltaSeed;
  class  HitData_t;

  namespace DeltaFinderTypes {

//-----------------------------------------------------------------------------
// diagnostics structure
//-----------------------------------------------------------------------------
    struct McPart_t : public TObject {
      const SimParticle*            fSim;             // type undefined here
      const DeltaCandidate*         fDelta;
      const ProtonCandidate*        fProton;
      std::vector<const HitData_t*> fListOfHits;
      int                           fFirstStation;
      int                           fLastStation;
      int                           fID;                // SimParticle::id().asInt()
      int                           fMotherID;
      int                           fPdgID;
      int                           fNHitsCE;           // N(hits) flagged as 'bkg'
      int                           fNChFlaggedDelta;   // (observable)
      int                           fNHitsDelta;        // N(hits)associated with reconstructed delta electron(s)

      int                           fNChFlaggedProton;  // number of hits flagged as proton (observable)
      int                           fNHitsProton;       // N(hits) associated with reconstructed proton(s)

      float                         fTime;
      float                         fHitTMin;           // min and max times of the straw hits
      float                         fHitTMax;

      float                         fStartMom;

      McPart_t(const SimParticle* Sim = nullptr): TObject() {
        fSim          = Sim;
        fDelta        = nullptr;
        fProton       = nullptr;
        fFirstStation = 999;
        fLastStation  = -1;
        fID           = -1;
        fMotherID     = -1;
        fPdgID            = 0;
        fNHitsCE          = 0;
        fNChFlaggedDelta  = 0;
        fNHitsDelta       = 0;
        fNChFlaggedProton = 0;
        fNHitsProton      = 0;
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
      int   PdgID    () { return fPdgID; }
    };

  }
}
#endif
