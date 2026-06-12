#ifndef CaloConditions_CaloSimCrystalConfig_hh
#define CaloConditions_CaloSimCrystalConfig_hh
//
// Initialize CaloMCCalib from fcl
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include <string>

namespace mu2e {

  struct CalSimCrystalConfig {
    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Atom<std::string> fileName{Name("fileName"),Comment("Calo MC calib filename")};
    fhicl::Atom<int>         verbose {Name("verbose"), Comment("verbosity")};
    fhicl::Atom<bool>        useDb   {Name("useDb"),   Comment("use database or fcl")};
  };

}

#endif
