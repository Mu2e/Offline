#ifndef TrackerConditions_FullReadoutStrawMaker_hh
#define TrackerConditions_FullReadoutStrawMaker_hh
//
// Initialize a list of straws with full readout
//

#include "Offline/TrackerConditions/inc/FullReadoutStraw.hh"
#include "Offline/DataProducts/inc/PanelId.hh"
#include "Offline/DataProducts/inc/PlaneId.hh"
#include "Offline/TrackerConfig/inc/FullReadoutStrawConfig.hh"
#include <set>
#include <iosfwd>

namespace mu2e {

  class FullReadoutStrawMaker {

    public:
      FullReadoutStrawMaker(FullReadoutStrawConfig const& config):_config(config) {}
      FullReadoutStraw::ptr_t fromFcl();
      FullReadoutStraw::ptr_t fromDb( /* db tables will go here*/ );

    private:

      // this object needs to be thread safe,
      // _config should only be initialized once
      const FullReadoutStrawConfig _config;
  };


} // namespace mu2e

#endif /* TrackerConditions_FullReadoutStrawMaker_hh */
