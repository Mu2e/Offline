#ifndef __TrkPatRec_RobustHelixFinder_types_hh__
#define __TrkPatRec_RobustHelixFinder_types_hh__

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

//#include "RecoDataProducts/inc/KalSeed.hh"

namespace mu2e {

  class RobustHelixFit;

  namespace RobustHelixFinderTypes {
  
    struct Data_t {
      const art::Event*               event;
      RobustHelixFit*                 result;
      fhicl::ParameterSet*            timeOffsets;

      enum  { kMaxHelicities = 2, kMaxSeeds = 100 };
      
      int     nTimePeaks;               // number of time peaks (input)
      int     nseeds   [kMaxHelicities]; // 0:all, 1:nhits > nhitsMin; assume nseeds <= 100
      int     ntclhits [kMaxHelicities][kMaxSeeds];
      int     nhits    [kMaxHelicities][kMaxSeeds];
      int     ntriplet0[kMaxHelicities][kMaxSeeds];
      int     ntriplet1[kMaxHelicities][kMaxSeeds];
      int     ntriplet2[kMaxHelicities][kMaxSeeds];
      
      double  r0       [kMaxHelicities][kMaxSeeds];
      double  radius   [kMaxHelicities][kMaxSeeds];
      double  lambda0  [kMaxHelicities][kMaxSeeds];
      double  lambda1  [kMaxHelicities][kMaxSeeds];
      double  chi2XY   [kMaxHelicities][kMaxSeeds];
      double  chi2ZPhi [kMaxHelicities][kMaxSeeds];
      double  pT       [kMaxHelicities][kMaxSeeds];
      double  p        [kMaxHelicities][kMaxSeeds];
      double  dr       [kMaxHelicities][kMaxSeeds];
      double  chi2d_helix[kMaxHelicities][kMaxSeeds];
      int maxSeeds() { return kMaxSeeds; }
      
    };
  }
}
#endif
