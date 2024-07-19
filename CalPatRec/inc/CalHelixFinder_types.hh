#ifndef CalPatRec_CalHelixFinder_types_hh
#define CalPatRec_CalHelixFinder_types_hh

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Sequence.h"

#include "art/Framework/Principal/Event.h"

#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTimeCluster.hh"
#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"
#include "Offline/Mu2eUtilities/inc/McUtilsToolBase.hh"

namespace art {
  class Event;
}

namespace fhicl {
  class ParameterSet;
}

class TH1F;
class TH2F;

namespace mu2e {

  namespace CalHelixFinderTypes {

    struct Config {
      fhicl::Atom<int>                      mcTruth   {fhicl::Name("mcTruth"),   fhicl::Comment("MC truth")                         };
      fhicl::Table<McUtilsToolBase::Config> mcUtils   {fhicl::Name("mcUtils"),   fhicl::Comment("MC Diag plugin")                   };
      fhicl::Atom<std::string>              tool_type {fhicl::Name("tool_type"), fhicl::Comment("tool type: Cal Helix Finder Diag") };
    };

    struct Data_t {
      const art::Event*    event;
      std::string          shLabel;
      //fhicl::ParameterSet* timeOffsets;

      enum  { kMaxSeeds = 100, kMaxHits = 200 };

      int     nTimePeaks;               // number of time peaks (input)
      int     nseeds   [        2]; // 0:all, 1:nhits > nhitsMin; assume nseeds <= 100
      int     ntclhits [kMaxSeeds];
      double  hitDzSeed[kMaxSeeds][kMaxHits];
      double  hitDrPred[kMaxSeeds][kMaxHits];
      int     nhits    [kMaxSeeds];
      int     good     [kMaxSeeds];
      double  radius   [kMaxSeeds];
      double  chi2XY   [kMaxSeeds];

      double  chi2ZPhi [kMaxSeeds];
      double  pT       [kMaxSeeds];
      double  p        [kMaxSeeds];
      int     nStationPairs[kMaxSeeds];
      double  dr       [kMaxSeeds];
      double  shmeanr  [kMaxSeeds];
      double  chi2d_helix[kMaxSeeds];
      double  npoints_loop0[kMaxSeeds];
      double  npoints_loop1[kMaxSeeds];
      double  chi2d_loop0[kMaxSeeds];
      double  chi2d_loop1[kMaxSeeds];
      double  chi2d_line_loop0[kMaxSeeds];
      double  chi2d_line_loop1[kMaxSeeds];
      int     loopId[kMaxSeeds];
      int maxSeeds() { return kMaxSeeds; }
      float   nHitsRatio[kMaxSeeds];
      float   eDepAvg[kMaxSeeds];
    };

    struct Hist_t {
      TH1F*  nTimePeaks;
      TH1F*  ntclhits[2];
      TH2F*  drVsDzSeed[2];
      TH1F*  nhits;           // number of hits on a helix
      TH1F*  nseeds  [2];
      TH1F*  radius  [2];
      TH1F*  chi2XY  [2];
      TH1F*  chi2ZPhi[2];
      TH1F*  pT      [2];
      TH1F*  p       [2];
      TH2F*  nhitsvspT;
      TH2F*  nhitsvsp;
      TH1F*  nStationPairs;
      TH1F*  dr[2];
      TH1F*  shmeanr[2];
      TH1F*  chi2d_helix[2];
      TH1F*  npoints_loop0;
      TH1F*  npoints_loop1;
      TH1F*  chi2d_loop0[2];
      TH1F*  chi2d_loop1[2];
      TH1F*  chi2d_line_loop0[2];
      TH1F*  chi2d_line_loop1[2];
      TH1F*  loopId[2];
      TH1F*  nHitsRatio[2];
      TH1F*  eDepAvg[2];
    };

  }

}
#endif
