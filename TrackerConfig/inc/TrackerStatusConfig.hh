#ifndef TrackerConditions_TrackerStatusConfig_hh
#define TrackerConditions_TrackerStatusConfig_hh
//
// Initialize TrackerStatus from fcl
//
#include <string>
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Tuple.h"

namespace mu2e {

  struct TrackerStatusSettings {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Atom<int> verbose{ Name("verbose"), Comment("verbosity: 0,1,2")};
    fhicl::Atom<bool> useDb{ Name("useDb"), Comment("use database or fcl")}; 
  };
  
  struct TrackerStatusConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    using TSSTable = fhicl::Table<TrackerStatusSettings>;
    using TSSequence = fhicl::Sequence<fhicl::Tuple<std::string,std::string,std::string>>;
    TSSTable settings { Name("Settings") };
    TSSequence status { Name("Status"),
      Comment("Provide in order: StrawId contained by the element in'plane_panel_straw' format \n"
	  "Level of the elemtn ('station', 'plane', 'panel', 'uniquepanel', 'straw', 'uniquestraw') \n"
	  "Colon-separated status bits to set for the element(see StrawStatus for details)") };
  };

}

#endif
