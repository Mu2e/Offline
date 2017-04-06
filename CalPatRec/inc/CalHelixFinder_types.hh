#ifndef __CalPatRec_CalHelixFinder_types_hh__
#define __CalPatRec_CalHelixFinder_types_hh__

#include "TObject.h"
#include "TH1.h"
#include "TH2.h"

namespace art {
  class Event;
};

namespace fhicl {
  class ParameterSet;
};

namespace mu2e {

  class CalHelixFinderData;

  struct CalHelixFinder_Data_t : public TObject {
    const art::Event*    event;
    CalHelixFinderData*  result;
    std::string          shLabel;
    fhicl::ParameterSet* timeOffsets;

    enum  { kMaxSeeds = 100 };
    int     nseeds  [2]; // 0:all, 1:nhits > nhitsMin; assume nseeds <= 100
    int     nhits   [kMaxSeeds];
    int     good    [kMaxSeeds];
    double  radius  [kMaxSeeds];
    double  chi2XY  [kMaxSeeds];
    double  chi2ZPhi[kMaxSeeds];
    double  pT      [kMaxSeeds];
    double  p       [kMaxSeeds];

    int maxSeeds() { return kMaxSeeds; }
  };

  struct CalHelixFinder_Hist_t : public TObject {
    TH1F*  nhits;           // number of hits on a helix  
    TH1F*  nseeds  [2];
    TH1F*  radius  [2];   
    TH1F*  chi2XY  [2];
    TH1F*  chi2ZPhi[2];
    TH1F*  pT      [2];
    TH1F*  p       [2];
    TH2F*  nhitsvspT;
    TH2F*  nhitsvsp;
  };

}
#endif
