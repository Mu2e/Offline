/*
#ifndef TrkPatRec_KalSeedFit_types_hh
#define TrkPatRec_KalSeedFit_types_hh

#include "TObject.h"
#include <vector>

class TH1F;
class TH2F;

namespace art {
  class Event;
}

namespace fhicl {
  class ParameterSet;
}

#include "Offline/RecoDataProducts/inc/KalSeed.hh"

namespace mu2e {

  class KalFitData;

  namespace KalSeedFitTypes {

    struct Data_t {
      const art::Event*               event;
      KalFitData*                     result;
      fhicl::ParameterSet*            timeOffsets;

      KalSeedCollection*              tracks;        // these report their momentum very unwillingly
      std::vector<int>                nrescued;      // per track
      std::vector<float>              mom;
    };
  }
}
#endif
*/
