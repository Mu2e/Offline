#ifndef TrackerConditions_Mu2eDetectorConfig_hh
#define TrackerConditions_Mu2eDetectorConfig_hh
//
// Initialize Mu2eDetector from fcl
//
#include <string>
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

namespace mu2e {

  struct Mu2eDetectorConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Atom<int> verbose{
      Name("verbose"), Comment("verbosity: 0,1,2")};
    fhicl::Atom<bool> useDb{
      Name("useDb"), Comment("use database or fcl")}; 
  };

}

#endif
