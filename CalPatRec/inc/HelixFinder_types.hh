#ifndef CalPatRec_HelixFinder_types_hh
#define CalPatRec_HelixFinder_types_hh

#include "fhiclcpp/types/Atom.h"

namespace mu2e {

namespace HelixFinderTypes {

struct Config {
  fhicl::Atom<std::string> tool_type{fhicl::Name("tool_type"),
                                     fhicl::Comment("tool type: HelixFinderDiag")};
};

struct tcInfo {
  int nHelices;
  int nComboHits;
  int nStrawHits;
  float time;
};

struct diagInfo {
  float moduleTime;
  int nHelices;
  int nTimeClusters;
  std::vector<tcInfo> timeClusterData;
};

} // namespace HelixFinderTypes
} // namespace mu2e
#endif
