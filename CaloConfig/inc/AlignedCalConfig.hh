#ifndef CaloConfig_AlignedCalConfig_hh
#define CaloConfig_AlignedCalConfig_hh
//
// Initialize AlignedCal from fcl
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include <string>

namespace mu2e {

  struct AlignedCalConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Atom<std::string> filenameDisk {
      Name("filenameDisk"), Comment("Filename for disk alignment file")};
    fhicl::Atom<int> verbose{
      Name("verbose"), Comment("verbosity")};
    fhicl::Atom<bool> useDb{
      Name("useDb"), Comment("use database or fcl")};
  };
}

#endif
