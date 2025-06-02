#ifndef CalPatRec_AgnosticHelixFinder_types_hh
#define CalPatRec_AgnosticHelixFinder_types_hh

#include "fhiclcpp/types/Atom.h"

namespace mu2e {

  namespace AgnosticHelixFinderTypes {

    struct Config {
      fhicl::Atom<std::string> tool_type{fhicl::Name("tool_type"), fhicl::Comment("tool type: AgnosticHelixFinderDiag")};
    };

    struct hitInfo {
      float eDep;
    };


    struct tcInfo {
      int     nHelices;
      int     nComboHits;
      int     nStrawHits;
    };

    struct hsInfo {
      float eDepAvg;
    };

    struct lineSegmentInfo {
      float   chi2dof;
      float   maxHitGap;
    };

    struct finalLineInfo{
      float nHitsRatio;
    };

    struct diagInfo {
      int                            nHelices;
      int                            nTimeClusters;
      std::vector<tcInfo>            timeClusterData;
      std::vector<hsInfo>            helixSeedData;
      std::vector<lineSegmentInfo>   lineSegmentData;
      std::vector<finalLineInfo>     lineInfoData;
      std::vector<hitInfo>           allHitsData;
    };

  } // namespace AgnosticHelixFinderTypes
} // namespace mu2e
#endif
