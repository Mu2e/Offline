#ifndef TrackerConditions_TrackerPanelMapConfig_hh
#define TrackerConditions_TrackerPanelMapConfig_hh
// clang-format off
// Initialize TrackerPanelMap from fcl
//
#include <string>
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Tuple.h"

namespace mu2e {

  struct TrackerPanelMapConfig {
    using Name   =fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Atom<int>     verbose{Name("verbose"), Comment("verbosity: 0,1,2")};
    fhicl::Atom<bool>    useDb  {Name("useDb"  ), Comment("use database or fcl")};

    fhicl::Sequence<int> mnid   {Name("mnid"   ), Comment("minnesota ID of the panel, main index")};
    fhicl::Sequence<int> dtcid  {Name("dtcid"  ), Comment("ID of the DTC reading the panel")};
    fhicl::Sequence<int> link   {Name("link"   ), Comment("DTC link of the panel (0-5)")};
    fhicl::Sequence<int> uniquePlane  {Name("uniquePlane"  ), Comment("geo index of the plane for this panel")};
    fhicl::Sequence<int> ppid   {Name("ppid"   ), Comment("production index of the plane")};
    fhicl::Sequence<int> panel  {Name("panel"  ), Comment("geo index of the panel")};
    fhicl::Sequence<int> zface  {Name("zface"  ), Comment("z-ordered index of this panel\'s face (0-3)")};
  };
}

#endif
