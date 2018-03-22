#ifndef __CalPatRec_CalSeedFit_types_hh__
#define __CalPatRec_CalSeedFit_types_hh__

#include "TObject.h"
#include <vector>

class TH1F;
class TH2F;

namespace art {
  class Event;
};

namespace fhicl {
  class ParameterSet;
};

#include "RecoDataProducts/inc/KalSeed.hh"

namespace mu2e {

  class KalFitResultNew;

  namespace CalSeedFitTypes {
  
    struct Data_t {
      const art::Event*               event;
      KalFitResultNew*                result;
      fhicl::ParameterSet*            timeOffsets;
      
      KalSeedCollection*              tracks;        // these report their momentum very unwillingly
      std::vector<int>                nrescued;	     // per track
      std::vector<float>              mom;
    };
  }
}
#endif
