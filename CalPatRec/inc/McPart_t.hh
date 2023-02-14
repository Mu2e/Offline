///////////////////////////////////////////////////////////////////////////////
// ComboHit.hh needs the definition of ProductID ...
///////////////////////////////////////////////////////////////////////////////
#ifndef __CalPatRec_inc_MCPart_t_hh__
#define __CalPatRec_inc_MCPart_t_hh__

#include "TObject.h"

namespace mu2e {

  class  Panel;
  class  SimParticle;
  struct DeltaCandidate;
  class  DeltaSeed;
  class  HitData_t;

  namespace DeltaFinderTypes {

//-----------------------------------------------------------------------------
// diagnostics structure
//-----------------------------------------------------------------------------
    struct McPart_t : public TObject {
      const SimParticle*            fSim;        // type undefined here
      const DeltaCandidate*         fDelta;
      std::vector<const HitData_t*> fListOfHits;
      int                           fFirstStation;
      int                           fLastStation;
      int                           fID;         // SimParticle::id().asInt()
      int                           fMotherID;
      int                           fPdgID;
      int                           fNHitsCE;  // N(hits) flagged as 'bkg'
      int                           fNHitsFlaggedBkg;
      int                           fNHitsDelta; // number of hits associated with reconstructed delta electrons
      float                         fTime;
      float                         fHitTMin;    // min and max times of the straw hits
      float                         fHitTMax;
      float                         fStartMom;

      McPart_t(const SimParticle* Sim = NULL): TObject() {
        fSim          = Sim;
        fDelta        = NULL;
        fID           = -1;
        fMotherID     = -1;
        fPdgID        = 0;
        fNHitsDelta   = 0;
        fNHitsCE      = 0;
        fNHitsFlaggedBkg = 0;
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

  };
}
#endif
