///////////////////////////////////////////////////////////////////////////////
// ComboHit.hh needs the definition of ProductID ...
///////////////////////////////////////////////////////////////////////////////
#ifndef CalPatRec_DeltaFinder_structures_hh
#define CalPatRec_DeltaFinder_structures_hh

#include "Offline/CalPatRec/inc/CalPatRec_enums.hh"
#include "Offline/CalPatRec/inc/ChannelID.hh"
#include "Offline/CalPatRec/inc/HitData_t.hh"

namespace mu2e {

  class  Panel;
  class  SimParticle;
  struct DeltaCandidate;
  class  DeltaSeed;

  namespace DeltaFinderTypes {

    extern float  stationZ   [kNStations];

    struct PhiPrediction_t {
      float fPhi;                       // predicted phi itself, -100 if no prediction
      float fErr;                       // uncertainty, defines the window
      int   fPanelID;                   // if can predict the panel, -1 otherwise
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
  }
}
#endif
