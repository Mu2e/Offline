#ifndef SimulationConfig_SimBookkeeperConfig_hh
#define SimulationConfig_SimBookkeeperConfig_hh

//
// Configuration for the MVACatalog ProditionsEntitiy
//

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Sequence.h"

namespace mu2e {

  struct SimStageEffConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Atom<std::string> tag{
      Name("tag"), Comment("Identifying tag for this efficiency")};
    fhicl::Atom<double> eff{
      Name("eff"), Comment("EFficiency value")};
  };

  struct SimBookkeeperConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Atom<int> verbose{
      Name("verbose"), Comment("verbosity: 0 or 1")}; 
    fhicl::Atom<bool> useDb{
      Name("useDb"), Comment("use database or fcl")}; 
    fhicl::Sequence< fhicl::Table<mu2e::SimStageEffConfig>> simStageEfficiencies{
      Name("simStageEfficiencies"), Comment("Sequence of SimStageEfficiency configurations")};
  };
}

#endif
