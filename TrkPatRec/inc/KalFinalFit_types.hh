#ifndef __TrkPatRec_KalFinalFit_types_hh__
#define __TrkPatRec_KalFinalFit_types_hh__

#include "TObject.h"
#include "TH1.h"

#include <vector>
#include "RecoDataProducts/inc/Doublet.hh"
#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/KalSeed.hh"

namespace art {
  class Event;
};

namespace mu2e {
  
  class KalFitData;
  class TTracker;
  class Calorimeter;
  class DoubletAmbigResolver;

  namespace KalFinalFitTypes {
  
    struct Data_t {
      const art::Event*     event;
      const TTracker*       tracker;
      const Calorimeter*    calorimeter;
      fhicl::ParameterSet*  timeOffsets;

      KalFitData*           result;
      KalRepCollection*     tracks;
      KalSeedCollection*    kscol;
      std::vector<Doublet>* listOfDoublets;
      DoubletAmbigResolver* dar;
      int                   eventNumber;
    };
  }
}
#endif
