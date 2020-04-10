#ifndef AnalysisConfig_MVACatalogConfig_hh
#define AnalysisConfig_MVACatalogConfig_hh

//
// Configuration for the MVACatalog ProditionsEntitiy
//

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Sequence.h"

namespace mu2e {
  struct MVAEntryConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;

    fhicl::Atom<std::string> trainName{
      Name("trainName"), Comment("Name of this training")};
    fhicl::Atom<std::string> xmlFileName{
      Name("xmlFileName"), Comment("XML file name for this training")};
    fhicl::Atom<bool> calibrated{
      Name("calibrated"), Comment("MVA training calibrated?"), false};
  };

  struct MVACatalogConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Atom<int> verbose{
      Name("verbose"), Comment("verbosity: 0 or 1")}; 
    fhicl::Atom<bool> useDb{
      Name("useDb"), Comment("use database or fcl")}; 
    fhicl::Sequence< fhicl::Table<MVAEntryConfig> > mvaConfigs{
      Name("mvaConfigs"), Comment("list of configurations of available MVA trainings")};
  };
}

#endif
