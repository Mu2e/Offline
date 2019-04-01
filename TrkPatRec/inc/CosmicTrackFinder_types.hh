#ifndef __TrkPatRec_CosmicTrackFinder_types_hh__
#define __TrkPatRec_CosmicTrackFinder_types_hh__

#include "TObject.h"
#include <vector>

namespace art {
  class Event;
};

namespace fhicl {
  class ParameterSet;
};

namespace mu2e {

  class CosmicTrackFit;

  namespace CosmicTrackFinderTypes{
  
    struct Data_t {
      const art::Event*               event;
      fhicl::ParameterSet*            timeOffsets;
      enum {kMaxSeeds = 1000, kMaxNHits = 10 };
      int     nseeds ; 
      int     npasses ;
      int     chi2numbers ;
      int     nTimePeaks;               // number of time peaks (input)
      
      int     nChPPanel[kMaxSeeds];    // maximum number of combohits per panel found in the TimeCluster
      int     nChHits[kMaxSeeds];    // maximum number of combohits per panel found in the TimeCluster
      int     ntclhits [kMaxSeeds]; //number of time cluster hits for each seed
      size_t  nhits    [kMaxSeeds]; //number of hits for each seed
      int     ninters  [kMaxSeeds];  
      double  chi2d_track[kMaxSeeds]; //chi-squared of fit for each tracl seed
      int     niters    [kMaxSeeds];
      int     nShFit [kMaxSeeds];
      int     nChFit [kMaxSeeds];
      int     good     [kMaxSeeds];
      int     nXYSh    [kMaxSeeds];
      int     nXYCh    [kMaxSeeds];
      //for chi2 fitting:
      double  hit_residualX[kMaxSeeds][kMaxNHits];
      double  hit_errorX [kMaxSeeds][kMaxNHits];
      double  hit_residualY[kMaxSeeds][kMaxNHits];
      double  hit_errorY [kMaxSeeds][kMaxNHits];
        
      int     nsh    [kMaxSeeds];
      double  chi2  [kMaxSeeds];
            
      int maxSeeds() { return kMaxSeeds; }
    };
   
  }
}
#endif
