#ifndef CaloConditions_CaloDAQMapConfig_hh
#define CaloConditions_CaloDAQMapConfig_hh
//
// Initialize CaloDAQMap from fcl
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include <string>

namespace mu2e {

  struct CaloDAQMapConfig {
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
