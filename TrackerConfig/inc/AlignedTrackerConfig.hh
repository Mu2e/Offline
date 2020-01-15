#ifndef TrackerConditions_AlignedTrackerConfig_hh
#define TrackerConditions_AlignedTrackerConfig_hh
//
// Initialize AlignedTracker from fcl
//
#include <string>
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

namespace mu2e {

  struct AlignedTrackerConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Atom<int> verbose{
      Name("verbose"), Comment("verbosity: 0,1,2")};
    fhicl::Atom<bool> useDb{
      Name("useDb"), Comment("use database or fcl")}; 
  };

}

#endif
