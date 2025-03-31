#ifndef TrkPatRec_RobustHelixFinder_types_hh
#define TrkPatRec_RobustHelixFinder_types_hh

#include "TObject.h"
#include <vector>
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

class TH1F;
class TH2F;

namespace art {
  class Event;
}

namespace fhicl {
  class ParameterSet;
}

//#include "RecoDataProducts/inc/KalSeed.hh"

namespace mu2e {

  class RobustHelixFit;

  namespace RobustHelixFinderTypes {

    struct Config {
      fhicl::Atom<int>         diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic level")};
      fhicl::Atom<std::string> tool_type{fhicl::Name("tool_type"), fhicl::Comment("tool type: Robust Helix Finder Diag")};
      fhicl::Atom<int>         mcTruth{fhicl::Name("mcTruth"), fhicl::Comment("MC truth")};
    };

    struct Data_t {
      const art::Event*               event;
      //      RobustHelixFit*                 result;
      fhicl::ParameterSet*            timeOffsets;

      enum  { kMaxHelicities = 2, kMaxSeeds = 100, kMaxNHits = 500 };

      int     nTimePeaks;               // number of time peaks (input)
      int     nChPPanel[kMaxHelicities][kMaxSeeds];    // maximum number of combohits per panel found in the TimeCluster
      int     nChHits  [kMaxHelicities][kMaxSeeds];    // maximum number of combohits per panel found in the TimeCluster
      int     nseeds   [kMaxHelicities]; // 0:all, 1:nhits > nhitsMin; assume nseeds <= 100
      int     ntclhits [kMaxHelicities][kMaxSeeds];
      int     nhits    [kMaxHelicities][kMaxSeeds];
      int     ntriplet0[kMaxHelicities][kMaxSeeds];
      int     ntriplet1[kMaxHelicities][kMaxSeeds];
      int     ntriplet2[kMaxHelicities][kMaxSeeds];
      int     xyniter  [kMaxHelicities][kMaxSeeds];
      int     fzniter  [kMaxHelicities][kMaxSeeds];
      int     niter    [kMaxHelicities][kMaxSeeds];
      int     nShFitXY [kMaxHelicities][kMaxSeeds];
      int     nChFitXY [kMaxHelicities][kMaxSeeds];
      int     nShFitCircle [kMaxHelicities][kMaxSeeds];
      int     nChFitCircle [kMaxHelicities][kMaxSeeds];
      int     nrescuedhits[kMaxHelicities][kMaxSeeds];
      int     nXYSh    [kMaxHelicities][kMaxSeeds];
      int     nZPhiSh  [kMaxHelicities][kMaxSeeds];

      int     nshsxy_0    [kMaxHelicities][kMaxSeeds];
      double  rsxy_0      [kMaxHelicities][kMaxSeeds];
      double  chi2dsxy_0  [kMaxHelicities][kMaxSeeds];

      int     nshsxy_1    [kMaxHelicities][kMaxSeeds];
      double  rsxy_1      [kMaxHelicities][kMaxSeeds];
      double  chi2dsxy_1  [kMaxHelicities][kMaxSeeds];

      int     nfz0counter [kMaxHelicities][kMaxSeeds];

      int     nshszphi_0  [kMaxHelicities][kMaxSeeds];
      double  lambdaszphi_0[kMaxHelicities][kMaxSeeds];
      double  chi2dszphi_0[kMaxHelicities][kMaxSeeds];

      int     nshszphi_1  [kMaxHelicities][kMaxSeeds];
      double  lambdaszphi_1[kMaxHelicities][kMaxSeeds];
      double  chi2dszphi_1[kMaxHelicities][kMaxSeeds];

      double  rinit    [kMaxHelicities][kMaxSeeds];
      double  radius   [kMaxHelicities][kMaxSeeds];
      double  lambda0  [kMaxHelicities][kMaxSeeds];
      double  lambda1  [kMaxHelicities][kMaxSeeds];
      double  chi2XY   [kMaxHelicities][kMaxSeeds];
      double  chi2ZPhi [kMaxHelicities][kMaxSeeds];
      double  pT       [kMaxHelicities][kMaxSeeds];
      double  p        [kMaxHelicities][kMaxSeeds];
      double  dr       [kMaxHelicities][kMaxSeeds];
      double  chi2d_helix[kMaxHelicities][kMaxSeeds];

      double  hitDr    [kMaxHelicities][kMaxSeeds][kMaxNHits];
      double  hitRWDot [kMaxHelicities][kMaxSeeds][kMaxNHits];

      int     nXYCh    [kMaxHelicities][kMaxSeeds];
      int     nZPhiCh  [kMaxHelicities][kMaxSeeds];

      int     nLoops   [kMaxHelicities][kMaxSeeds];

      float   meanHitRadialDist [kMaxHelicities][kMaxSeeds];

      int     nHitsLoopFailed [kMaxHelicities][kMaxSeeds];

      float   eDepAvg[kMaxHelicities][kMaxSeeds];

      int maxSeeds() { return kMaxSeeds; }

    };
  }
}
#endif
