#ifndef TrackerConditions_DeadStrawMaker_hh
#define TrackerConditions_DeadStrawMaker_hh
//
// Initialize a list of dead straws
//

// Mu2e includes
#include "TrackerConditions/inc/DeadStraw.hh"
#include "DataProducts/inc/PanelId.hh"
#include "DataProducts/inc/PlaneId.hh"
#include "TrackerConfig/inc/DeadStrawConfig.hh"

// C++ includes
#include <set>
#include <iosfwd>

namespace mu2e {

  class DeadStrawMaker {

  public:
    DeadStrawMaker(DeadStrawConfig const& config):_config(config) {}
    DeadStraw::ptr_t fromFcl();
    DeadStraw::ptr_t fromDb( /* db tables will go here*/ );

  private:

    void addDeadPlane( PlaneId const& id, DeadStraw::set_t& dead);
    void addDeadPanel( PanelId const& id, DeadStraw::set_t& dead);
    void addDeadStraw( StrawId const& id, DeadStraw::set_t& dead, 
		       double range=-1.0);

    // this object needs to be thread safe, 
    // _config should only be initialized once
    const DeadStrawConfig _config;
  };


} // namespace mu2e

#endif /* TrackerConditions_DeadStrawMaker_hh */
