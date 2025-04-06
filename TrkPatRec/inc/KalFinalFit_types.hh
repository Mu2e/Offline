/*
#ifndef TrkPatRec_KalFinalFit_types_hh
#define TrkPatRec_KalFinalFit_types_hh

#include "TObject.h"
#include "TH1.h"

#include <vector>
#include "Offline/BTrkData/inc/Doublet.hh"
#include "Offline/RecoDataProducts/inc/KalRepCollection.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"

namespace art {
  class Event;
}

namespace mu2e {

  class KalFitData;
  class Tracker;
  class Calorimeter;
  class DoubletAmbigResolver;

  namespace KalFinalFitTypes {

    struct Data_t {
      const art::Event*     event;
      const Tracker*       tracker;
      const Calorimeter*    calorimeter;
      fhicl::ParameterSet*  timeOffsets;

      KalFitData*           result;
      KalRepCollection*     tracks;
      KalSeedCollection*    kscol;
      std::vector<Doublet>* listOfDoublets;
      DoubletAmbigResolver* dar;
      int                   eventNumber;

      unsigned              tchDiskId;
      unsigned              tchAdded;
      double                tchDepth;
      double                tchDOCA;
      double                tchDt;
      double                tchTrkPath;
      double                tchEnergy;
    };
  }
}
#endif
*/
