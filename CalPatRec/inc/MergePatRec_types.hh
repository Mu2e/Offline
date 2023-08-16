#ifndef CalPatRec_MergePatRec_types_hh
#define CalPatRec_MergePatRec_types_hh

#include "Offline/RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"

namespace art {
  class Event;
}

namespace fhicl {
  class ParameterSet;
}

namespace mu2e {
  class Tracker;

  namespace MergePatRecTypes {
//-----------------------------------------------------------------------------
// data structure passed to the histogramming routine
//-----------------------------------------------------------------------------
    struct Data_t {
      const art::Event*       event;
      const Tracker*          tracker;
      const KalRepPtrCollection*    list_of_kreps_tpr;
      const KalRepPtrCollection*    list_of_kreps_cpr;
      const KalSeedCollection*      list_of_kseed_tpr;
      const KalSeedCollection*      list_of_kseed_cpr;
      int                     debugLevel;             // printout level
    };
  }
}
#endif
