#ifndef CalPatRec_AgnosticHelixFinder_types_hh
#define CalPatRec_AgnosticHelixFinder_types_hh

#include "fhiclcpp/types/Atom.h"

namespace mu2e {

  namespace AgnosticHelixFinderTypes {

    struct Config {
      fhicl::Atom<std::string> tool_type{fhicl::Name("tool_type"), fhicl::Comment("tool type: AgnosticHelixFinderDiag")};
    };

    struct tcInfo {
      int     nHelices;
      int     nComboHits;
      int     nStrawHits;
      float   time;
    };

    struct lineSegmentInfo {
      float   chi2dof;
      float   maxHitGap;
    };

    struct diagInfo {
      float                          moduleTime;
      int                            nHelices;
      int                            nTimeClusters;
      std::vector<tcInfo>            timeClusterData;
      std::vector<lineSegmentInfo>   lineSegmentData;
    };

  } // namespace AgnosticHelixFinderTypes
} // namespace mu2e
#endif
