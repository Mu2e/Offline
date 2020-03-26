#ifndef TrkDiag_TrkQualCatalogConfig_hh
#define TrkDiag_TrkQualCatalogConfig_hh

//
// Configuration for the TrkQualCatalog ProditionsEntitiy
//

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Sequence.h"

namespace mu2e {
  struct TrkQualEntryConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;

    fhicl::Atom<std::string> trainName{
      Name("trainName"), Comment("Name of this training")};
    fhicl::Atom<std::string> xmlFileName{
      Name("xmlFileName"), Comment("XML file name for this training")};
    fhicl::Atom<bool> calibrated{
      Name("calibrated"), Comment("TrkQual training calibrated?"), false};
  };

  struct TrkQualCatalogConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Atom<int> verbose{
      Name("verbose"), Comment("verbosity: 0 or 1")}; 
    fhicl::Atom<bool> useDb{
      Name("useDb"), Comment("use database or fcl")}; 
    fhicl::Sequence< fhicl::Table<TrkQualEntryConfig> > trkQualConfigs{
      Name("trkQualConfigs"), Comment("list of configurations of available TrkQual trainings")};
  };
}

#endif
