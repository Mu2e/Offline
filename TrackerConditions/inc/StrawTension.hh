#ifndef TrackerConditions_StrawTension_hh
#define TrackerConditions_StrawTension_hh
//
// Class to specify the initial wire and straw tensions.  Actual tensions are modeled as an exponential decay
// 
// Original author David Brown (7/2020)
//
// Mu2e includes
#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/StrawEnd.hh"
#include "canvas/Persistency/Provenance/Timestamp.h"
namespace mu2e {
  struct StrawTension {
    StrawId id_; // defines which straw.
    float tension_; // tension in gm;
    art::Timestamp measure_; // time at which the tension was measured
  };
}
#endif
