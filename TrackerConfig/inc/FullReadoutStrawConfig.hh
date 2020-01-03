#ifndef TrackerConditions_FullReadoutStrawConfig_hh
#define TrackerConditions_FullReadoutStrawConfig_hh
//
// Initialize FullReadoutStraw from fcl
//
#include <string>
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

namespace mu2e {

  struct FullReadoutStrawConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Atom<int> verbose{
      Name("verbose"), Comment("verbosity: 0,1,2")};
    fhicl::Atom<bool> useDb{
      Name("useDb"), Comment("use database or fcl")}; 
    fhicl::Sequence<std::string> straws{
      Name("straws"), Comment("strawId's as strings") };
  };

}

#endif
