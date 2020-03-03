#ifndef TrackerConditions_DeadStrawConfig_hh
#define TrackerConditions_DeadStrawConfig_hh
//
// Initialize DeadStraw from fcl
//
#include <string>
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

namespace mu2e {

  struct DeadStrawConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Atom<int> verbose{
      Name("verbose"), Comment("verbosity: 0,1,2")};
    fhicl::Atom<bool> useDb{
      Name("useDb"), Comment("use database or fcl")}; 
    fhicl::Sequence<int> deadPlanes{
      Name("deadPlanes"), Comment("int sequence") };
    fhicl::Sequence<std::string> deadPanels{
      Name("deadPanels"), Comment("string sequence") };
    fhicl::Sequence<std::string> deadStraws{
      Name("deadStraws"), Comment("string sequence") };
    fhicl::Sequence<std::string> partlyDeadStraws{
      Name("partlyDeadStraws"), Comment("string sequence") };
  };

}

#endif
