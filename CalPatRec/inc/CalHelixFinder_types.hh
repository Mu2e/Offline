#ifndef __CalPatRec_CalHelixFinder_types_hh__
#define __CalPatRec_CalHelixFinder_types_hh__

namespace art {
  class Event;
};

namespace fhicl {
  class ParameterSet;
};

class TH1F;
class TH2F;

namespace mu2e {

  namespace CalHelixFinderTypes {
    
    struct Data_t {
      const art::Event*    event;
      std::string          shLabel;
      fhicl::ParameterSet* timeOffsets;

      enum  { kMaxSeeds = 100 };

      int     nseeds  [        2]; // 0:all, 1:nhits > nhitsMin; assume nseeds <= 100
      int     nhits   [kMaxSeeds];
      int     good    [kMaxSeeds];
      double  radius  [kMaxSeeds];
      double  chi2XY  [kMaxSeeds];
      double  chi2ZPhi[kMaxSeeds];
      double  pT      [kMaxSeeds];
      double  p       [kMaxSeeds];

      int maxSeeds() { return kMaxSeeds; }
    };

    struct Hist_t {
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

}
#endif
