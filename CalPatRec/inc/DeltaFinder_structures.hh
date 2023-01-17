#ifndef __CalPatRec_DeltaFinder_structures_hh__
#define __CalPatRec_DeltaFinder_structures_hh__

#include "Offline/CalPatRec/inc/DeltaFinder_enums.hh"

#include "Offline/RecoDataProducts/inc/ComboHit.hh"

namespace mu2e {

  namespace DeltaFinderTypes {
    extern float  stationZ   [kNStations];
  };

  class DeltaSeed;
  class Panel;
  class SimParticle;

  struct ChannelID {
    int Station;
    int Plane;
    int Face;
    int Panel;
    int Layer;
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
  struct DeltaSeed;

  struct HitData_t {
    const ComboHit*         fHit;
    DeltaSeed*              fSeed;           // nullptr if not associated...
       //      int                     fNSecondHits;
    int                     fDeltaIndex;
    float                   fChi2Min;
    float                   fSigW2;    // cached resolution^2 along the wire
    float                   fCorrTime; // cache hit corrected time

    HitData_t(const ComboHit* Hit) {
      fHit         = Hit;
      fChi2Min     = 99999.0;
      float sigw   =  Hit->posRes(ComboHit::wire);
      fSigW2       = sigw*sigw;  // to be used to calculate chi2...
      fSeed        = nullptr;
      // fNSecondHits =  0;
      fDeltaIndex  = -1;
      fCorrTime    = Hit->correctedTime();
    }

    int Used() const { return (fSeed != nullptr) ; }
  };

  struct PanelZ_t {
    int                              fNHits  ; // guess, total number of ComboHits
    std::vector<HitData_t>*          fHitData;
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

}
#endif
