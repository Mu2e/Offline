#ifndef __CalPatRec_DeltaFinder_structures_hh__
#define __CalPatRec_DeltaFinder_structures_hh__

#include "TObject.h"

#include "Offline/CalPatRec/inc/DeltaFinder_enums.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"

namespace mu2e {

  class DeltaSeed;
  struct DeltaCandidate;

  class Panel;
  class SimParticle;

  namespace DeltaFinderTypes {
    extern float  stationZ   [kNStations];

    struct ChannelID {
      int Station;
      int Plane;
      int Face;
      int Panel;
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

    struct HitData_t {
      const ComboHit*         fHit;
      DeltaSeed*              fSeed;           // nullptr if not associated...
      int                     fUsed;           // TBD
      int                     fZFace;          // z-ordered face (for printing)
      float                   fChi2Min;
      float                   fSigW2;          // cached resolution^2 along the wire
      float                   fCorrTime;       // cached hit corrected time
      int                     fDeltaIndex;     // is it really needed? **FIXME**

      HitData_t(const ComboHit* Hit,int ZFace) {
        fHit         = Hit;
        fSeed        = nullptr;
        fUsed        = 0;
        fZFace       = ZFace;
        fChi2Min     = 99999.0;
        float sigw   =  Hit->posRes(ComboHit::wire);
        fSigW2       = sigw*sigw;
        fCorrTime    = Hit->correctedTime();
        fDeltaIndex  = -1;
      }

      int Used() const { return fUsed ; }
    };
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
      int                           fNHitsCE;
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
