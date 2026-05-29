#ifndef CaloConditions_CalDAQMapConfig_hh
#define CaloConditions_CalDAQMapConfig_hh
//
// Initialize CalDAQMap from fcl
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include <string>

namespace mu2e {

  struct CalDAQMapConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Atom<std::string> fileSpec{
      Name("fileSpec"), Comment("Filespec for DAQMap file")};
    fhicl::Atom<int> verbose{
      Name("verbose"), Comment("verbosity: 0,1,2")};
    fhicl::Atom<bool> useDb{
      Name("useDb"), Comment("use database or fcl")};

  };

}

#endif
