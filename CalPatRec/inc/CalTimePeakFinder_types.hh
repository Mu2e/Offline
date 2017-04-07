#ifndef __CalPatRec_CalTimePeakFinder_types_hh__
#define __CalPatRec_CalTimePeakFinder_types_hh__

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
  class CaloCluster;
  class TimeCluster;

  struct CalTimePeakFinder_Data_t : public TObject {
    enum  { kMaxSeeds = 100 };

    const art::Event*    event;
    // KalFitResultNew*     result;
    // fhicl::ParameterSet* timeOffsets;

    const CaloCluster*  cl;
    const TimeCluster*  timeCluster;

    int                 nseeds[2];
    int                 minNHits;

    int maxSeeds() { return kMaxSeeds; }
  };

  struct CalTimePeakFinder_Hist_t : public TObject {
    TH1F*  nhits;           // number of hits on a helix  
    TH1F*  energy[2];   
    TH1F*  time[2];
    TH2F*  time_vs_nhits;
    TH2F*  energy_vs_nhits;
    TH1F*  nseeds[2];
  };

}
#endif
