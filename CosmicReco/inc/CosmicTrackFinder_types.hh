#ifndef __TrkPatRec_CosmicTrackFinder_types_hh__
#define __TrkPatRec_CosmicTrackFinder_types_hh__
//Stores diagnostic information about cosmic seed fit (residuals, chi2 before and after first chi2 minimization application)
#include "TObject.h"
#include <vector>

namespace art {
  class Event;
  class Run;
};

namespace fhicl {
  class ParameterSet;
};

namespace mu2e {

  class CosmicTrackFit;

  namespace CosmicTrackFinderTypes{
    
    struct Data_t {
      const art::Run*               run;
      const art::Event*               event;
      
      enum {kMaxSeeds = 1000, kMaxNHits = 50 };
      int     nseeds ; //Number of Seeds
      int     npasses ; //Number passing all cuts in finder
      int     nTimePeaks;               // number of time clusters
      
      int     ntclhits ; //number of time cluster hits for each seed
      size_t  nhits    ; //number of hits for each seed
      
      double  Final_chi2d_track; //chi-squared of fit for each tracl seed
      double  Final_chi2dX_track; //chi-squared of fit for each tracl seed
      double  Final_chi2dY_track; //chi-squared of fit for each tracl seed
      double  Initial_chi2dX_track; //chi-squared of fit for each tracl seed
      double  Initial_chi2dY_track; //chi-squared of fit for each tracl seed
      double  Initial_chi2d_track; //chi-squared of fit for each tracl seed
      double  Change_chi2d_track; //total change
      double  Change_chi2d_trackX; 
      double  Change_chi2d_trackY; 
      double  True_chi2d_track;
      
      
      int     niters    ; //Number of iterations on seed
      int     nShFit ; //Number ofstraw hits in the fit
      int     nChFit ; //Number of ComboHits Fit
      int     good     ; //Number which are "good"
       
      //for chi2 fitting:
      double  Final_hit_residualX[kMaxNHits];
      double  Final_hit_errorX [kMaxNHits];
      double  Final_hit_residualY[kMaxNHits];
      double  Final_hit_errorY [kMaxNHits];
      double  Final_hit_pullX[kMaxNHits];
      double  Final_hit_pullY[kMaxNHits];
      //for chi2 fitting:
      double  Initial_hit_residualX[kMaxNHits];
      double  Initial_hit_errorX [kMaxNHits];
      double  Initial_hit_residualY[kMaxNHits];
      double  Initial_hit_errorY [kMaxNHits];
      double  Initial_hit_pullX[kMaxNHits];
      double  Initial_hit_pullY[kMaxNHits];
      //Errors:
      
      double    Final_ErrorsTot[kMaxNHits];
      double    Initial_ErrorsTot[kMaxNHits];
      //int     nsh    ;
      
      int maxSeeds() { return kMaxSeeds; }
    };
     
  }
}
#endif
