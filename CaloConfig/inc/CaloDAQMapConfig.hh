#ifndef CaloConditions_CaloDAQMapConfig_hh
#define CaloConditions_CaloDAQMapConfig_hh
//
// Initialize CaloDAQMap from fcl
//
#include <string>
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

namespace mu2e {

  struct CaloDAQMapConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Sequence<uint16_t> DIRAC2CaloMap{
      Name("DIRAC2CaloMap"), Comment("map from DIRAC channel to Calo roID")};
    fhicl::Sequence<uint16_t> Calo2DIRACMap{
      Name("Calo2DIRACMap"), Comment("map from Calo roID to DIRAC channel")};
    fhicl::Atom<int> verbose{
      Name("verbose"), Comment("verbosity: 0,1,2")};
    fhicl::Atom<bool> useDb{
      Name("useDb"), Comment("use database or fcl")};

  };

}

#endif
