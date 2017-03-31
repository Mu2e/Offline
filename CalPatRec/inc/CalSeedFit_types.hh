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

namespace mu2e {

  class KalFitResultNew;

  struct CalSeedFit_Data_t : public TObject {
    enum  { kMaxSeeds = 100 };

    const art::Event*    event;
    KalFitResultNew*     result;
    fhicl::ParameterSet* timeOffsets;

    int                 ntracks;           // 0:all, 1:nhits > nhitsMin; assume nseeds <= 100
    std::vector<int>    nrescued;

    int maxSeeds() { return kMaxSeeds; }
  };

  struct CalSeedFit_Hist_t : public TObject {
    TH1F*  nhits;           // number of hits on a htrack candidate
    TH1F*  seeddoca   [3];
    TH1F*  chi2       [2];
    TH1F*  p          [2];
    TH1F*  NpointsSeed[2]; //
    TH1F*  ntracks    [2];
  };

}
#endif
