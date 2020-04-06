#ifndef TrackerConditions_TrackerDAQConditionsConfig_hh
#define TrackerConditions_TrackerDAQConditionsConfig_hh
//
// Initialize TrackerDAQConditions from fcl
//
#include <string>
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

namespace mu2e {

  struct TrackerDAQConditionsConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Sequence<uint16_t> DRACStrawMap{
      Name("DRACStrawMap"), Comment("map from DRAC channel to straw number")};
    fhicl::Sequence<uint16_t> ROCPanelMap{
      Name("ROCPanelMap"), Comment("map from ROC panel to panel number")};
    fhicl::Atom<int> verbose{
      Name("verbose"), Comment("verbosity: 0,1,2")};
    fhicl::Atom<bool> useDb{
      Name("useDb"), Comment("use database or fcl")}; 

  };

}

#endif
