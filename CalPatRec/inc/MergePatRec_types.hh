#ifndef __CalPatRec_MergePatRec_types_hh__
#define __CalPatRec_MergePatRec_types_hh__

#include "RecoDataProducts/inc/KalRepPtrCollection.hh"

namespace art {
  class Event;
};

namespace fhicl {
  class ParameterSet;
};

namespace mu2e {
  class TTracker;

  namespace MergePatRecTypes {
//-----------------------------------------------------------------------------
// data structure passed to the histogramming routine
//-----------------------------------------------------------------------------
    struct Data_t {
      const art::Event*       event;
      const TTracker*         tracker;
      const KalRepPtrCollection*    list_of_kreps_tpr;
      const KalRepPtrCollection*    list_of_kreps_cpr;
      int                     debugLevel;	     // printout level
    };
  }
}
#endif
